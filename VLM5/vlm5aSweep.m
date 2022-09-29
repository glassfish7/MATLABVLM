function [CLpot,CDiPot,CL,CDi,CDiTreffz] = vlm5aSweep(sections,data,alphaIncRoot,alphaIncTip,phi,camberLine,bPanels,cPanels,alpha,beta,Uinf,rho,LE,wakeT,os)
%VLM2SweepAlpha Sweeps Alpha for a given planform with VLM2
%   data = [taper;rootC;AR;LambLE] column for each section


%Set number of panels and alpha
NLE = LE; %Number of rows of panels to apply suction analogy to starting with the Leading Edge


%Use Geometry Engine to Obtain Panels
[controlPointsY,controlPointsX,panelQuarterC,panelGeomY,panelGeomX,VcontrolPointsX,S,b,N,K] = geomEngine(sections,data,bPanels,cPanels);

realAR = b^2/S;
fprintf('The Aspect Ratio of this Wing is :%f\n',realAR);
fprintf('The Wing Area is: %f\n',S);


%Set Twist
twist = spline([panelGeomY(1) panelGeomY(length(panelGeomY)/2 +0.5) panelGeomY(end)],[alphaIncTip alphaIncRoot alphaIncTip],controlPointsY);
twist = repmat(twist,cPanels,1);
twist = twist(:);


CLpot = zeros(1,length(alpha));
CDiPot = zeros(1,length(alpha));
CDiTreffz = zeros(1,length(alpha));
CDi = zeros(1,length(alpha));
CL = zeros(1,length(alpha));


%Set Airfoil sections
%Set normalized control point locations
normControlPoints = linspace(0,1,cPanels+1);
normControlPoints  = normControlPoints + (3/(4*cPanels));
normControlPoints  = normControlPoints(1:end-1);

%Assign camber line to each chordwise wing section
%Take derivative with central difference
%3 approaches: 4-digit NACA, custom camber line, and no camber
if isa(camberLine,'char')
    if length(camberLine) ~= 4
        error('Please enter a 4 digit NACA airfoil');
    else
        camber = zeros(1,length(normControlPoints));
        camberD = zeros(1,length(normControlPoints));
        m = str2double(camberLine(1))/100;
        p = str2double(camberLine(2))/10;
        for i = 1:length(normControlPoints)
            if normControlPoints <= p
                camber(i) = (m/p^2)*(2*p*normControlPoints(i) - normControlPoints(i).^2);
                camberD(i) = (2*m/p^2)*(p - normControlPoints(i));
            else
                camber(i) = (m/(1-p)^2)*(1-2*p+2*p*normControlPoints(i) - normControlPoints(i).^2);
                camberD(i) = (2*m/(1-p)^2)*(p - normControlPoints(i));
            end
        end
    end
    camberMat  = repmat(camberD',1,2*bPanels);
    camberDrep = repmat(camberD',2*bPanels,1);
elseif camberLine == 0
    camber = 0;
    camberDrep = zeros(N,1);
    camberMat = zeros(cPanels,2*bPanels);
else
    
    %Custom camber line functionality. 
    %Multi section airfoil functionality to be added later

    [camberRows,camberCol] = size(camberLine);

    if camberRows ~= 2
        error('Missing coordinates. Cannot interpolate camber line')
    end

    if camberLine(1,end) > 1
        error('Normalize camber line coordinates')
    end

    camberD = zeros(1,length(normControlPoints));
    camber = spline(camberLine(1,:),camberLine(2,:),normControlPoints);

    %Prep Central Diff 
    centDiffCont = zeros(1,2*length(normControlPoints));
    count = 1;
    for i = 1:2:2*length(normControlPoints)
        centDiffCont(i) = normControlPoints(count)-0.01;
        centDiffCont(i+1) = normControlPoints(count)+0.01;
        count = count +1;
    end
    camberCentCont = spline(camberLine(1,:),camberLine(2,:),centDiffCont);
    
    %Take central diff
    count = 1;
    for i = 1:2:2*length(normControlPoints)
            camberD(count) = (camberCentCont(i+1)-camberCentCont(i))/0.02;
            count = count + 1;
    end
    
    camberDrep = repmat(camberD',2*bPanels,1);
    camberMat  = repmat(camberD',1,2*bPanels);
end



for d2 = 1:length(alpha)
%Generate Vorticity Matrix

%Panels are numbered sequentially from 1 to the total number of panels
%Each panel is numbered in Figure 2

%The row in the matrix corresponding to a particular panel contains all of
%the downwash contribution on that panel from every other panel in the
%system

%The column in the matrix corresponding to a particular panel contains all
%of the contribution of that particular panel to every other panel in the
%system


vorticityMatrix = zeros(N,N);

%Set counters to 1
panelCount =1;
contCount =1;


%Loop through each panel on the wing. Calculate downwash contributions ONTO this
%panel 
for i = 1:length(controlPointsY)
    for j = 1:cPanels
        
    
        %Loop through each panel on the wing. Calculate downwash contributions BY this
        %panel 
        for k = 1:length(controlPointsY)
            for m = 1:cPanels
                    %Bound Vortex: Vortex 1
                    [u1,v1,w1] = VORTEX(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),0,-1*(panelQuarterC(m,k)-controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                    
                    %Trailing vorticies into the wake
                    %Vortex 2: TE to wake
                    %Vortex 3: wake to TE
                    [u2,v2,w2] = VORTEX(wakeT*b - (-1*(panelGeomX(end,k+1)-os*b)),tan(beta)*(wakeT*b - (-1*(panelGeomX(end,k+1)-os*b))),tan(alpha(d2))*(wakeT*b - (-1*(panelGeomX(end,k+1)-os*b))),...
                        -1*(panelGeomX(end,k+1)-os*b-controlPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i),0);
                    
                    [u3,v3,w3] = VORTEX((-1*(panelGeomX(end,k)-os*b))- wakeT*b ,tan(beta)*((-1*(panelGeomX(end,k)-os*b))- wakeT*b),tan(alpha(d2))*((-1*(panelGeomX(end,k)-os*b))- wakeT*b),...
                        wakeT*b-(-1*controlPointsX(j,i)),panelGeomY(k)-tan(beta)*((-1*(panelGeomX(end,k)-os*b))- wakeT*b)-controlPointsY(i),tan(alpha(d2))*(wakeT*b-(-1*panelGeomX(end,k))));

                    %Trailing Vorticies on the wing
                    %Vortex 4: quarter chord to TE
                    %Vortex 5: TE to quarter chord
                    [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)-os*b) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-controlPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i),0);
                    [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k)-os*b)),0,0,(-1*(panelGeomX(end,k)-os*b))-(-1*controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                    
                    %Sum contributions from each segment in the u,v,w
                    %directions
                    u = u1 + u2 + u3 + u4 + u5;
                    v = v1 + v2 + v3 + v4 + v5; 
                    w = w1 + w2 + w3 + w4 + w5;
                    
                    %Apply effects of camber and dihedral to contributions
                    %before assigning to vorticity matrix. Note
                    %asymmetircal dihedral angle on each wing
                   if panelCount > N/2
                    vorticityMatrix(panelCount,contCount) = w*cos(phi)*cos(atan(camberMat(j,i))) - u*sin(atan(camberMat(j,i)))*cos(phi) - v*sin(phi);
                   else
                    vorticityMatrix(panelCount,contCount) = w*cos(phi)*cos(atan(camberMat(j,i))) - u*sin(atan(camberMat(j,i)))*cos(phi) - v*sin(-phi);
                   end
                   contCount = contCount + 1;    
            end
        end
        
        %Update counters
        panelCount = panelCount + 1;
        contCount =1;
    end
end



%Generate RHS of matrix equation with no penetration boundary condition.
%Includes effects from alpha, beta, phi, and camber
%RHS1 is left wing and RHS2 is right wing
%Concatinate vectors 
RHS1 = Uinf*(sin(atan(camberDrep(1:N/2,:)))*cos(phi)*cos(alpha(d2))*cos(beta)+...
    sin(beta)*sin(-phi)...
    -sin(alpha(d2))*cos(beta)*cos(atan(camberDrep(1:N/2,:)))*cos(phi)).*ones(N/2,1);
RHS2 = Uinf*(sin(atan(camberDrep(N/2+1:end,:)))*cos(phi)*cos(alpha(d2))*cos(beta)+...
    sin(beta)*sin(phi)-...
    sin(alpha(d2))*cos(beta)*cos(atan(camberDrep(N/2+1:end,:)))*cos(phi)).*ones(N/2,1);
RHS = [RHS1;RHS2];
%No dihedral version for debug
% RHS = (Uinf*(-sin(alpha+twist)+camberDrep)).*ones(N,1);

%Solve for vorticity
gamma = linsolve(vorticityMatrix,RHS);

%Place gamma vector back into array for easy interpretation of results
gammaMatrix = zeros(cPanels,2*sum(bPanels));
gammaCount =1;
for p = 1:2*sum(bPanels)
    for o = 1:cPanels
        gammaMatrix(o,p) = gamma(gammaCount);
        gammaCount = gammaCount + 1;
    end
end


%This section of code breaks the gammaMatrix back into sections 
%by identifying start and ends points of each section in the matrix
sPoints = zeros(2,sections);
ePoints = zeros(2,sections);
for i = 1:sections
    sPoints(1,i) = 1 + sum(bPanels(sections:-1:i))-bPanels(i);
    sPoints(2,i) = sPoints(1,i) + sum(bPanels(1:i))+sum(bPanels(1:i-1));

    ePoints(1,i) = sPoints(1,i) + bPanels(i)-1;
    ePoints(2,i) = sPoints(2,i) + bPanels(i)-1;
end
%Use the start and end points to apply Kutta-Joukowski to each section
%with the appropriate width(K) of each section
for t = 1:sections
          potentialLiftDist(:,[sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])  = ...
              rho*Uinf*gammaMatrix(:,[sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])*K(t);

end

% In order to apply the Kutta-Joukowski the bound vortices are projected
% normal to the freestream. NOTE: Sweep is already accounted for in the panel width K from
% geomEngine.m
for cPan = 1:cPanels
    potentialLiftDist(cPan,1:bPanels) = potentialLiftDist(cPan,1:bPanels)*cos(-beta);
    potentialLiftDist(cPan,bPanels+1:end) = potentialLiftDist(cPan,bPanels+1:end)*cos(beta);
end

%Calculate the induced drag produced by the wing in potential flow with Treffz Plane Analysis
gammaTreffz = sum(gammaMatrix);
wiTreffz = zeros(1,length(controlPointsY));
%Find downwash at the wake
for i = 1:length(controlPointsY)
                    wiTreffz(i) = sum(gammaTreffz.*(-1./(panelGeomY(1:end-1)-controlPointsY(i)) + 1./(panelGeomY(2:end)-controlPointsY(i))))/(4*pi);   
end
inducedDragTreffz = (rho)*sum(gammaTreffz.*wiTreffz);
CDiTreffz(d2) = 2*inducedDragTreffz/(rho*Uinf^2*S);

%Vortex Lift Section
%Calculate Zero Leading edge suction lift 
%Find Normal Force on each panel for zero-suction lift
%Also calculate the potential induced drag
normalForce = zeros(cPanels,2*sum(bPanels));
potIndDrag =  zeros(cPanels,2*sum(bPanels));
suctionNormalForce = zeros(cPanels,2*sum(bPanels));
panelCount = 1;

for i = 1:length(controlPointsY)
    for j = 1:cPanels
        wi = 0;
        %Calculate downwash on panel      
        for k = 1:length(controlPointsY)
            for m = 1:cPanels         
                    %Bound Vortex: Vortex 1
                    [u1,v1,w1] = VORTEX(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),0,-1*(panelQuarterC(m,k)-VcontrolPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                    
                    %Trailing vorticies into the wake
                    %Vortex 2: TE to wake
                    %Vortex 3: wake to TE
                    [u2,v2,w2] = VORTEX(wakeT*b - (-1*(panelGeomX(end,k+1)-os*b)),tan(beta)*(wakeT*b - (-1*(panelGeomX(end,k+1)-os*b))),tan(alpha(d2))*(wakeT*b - (-1*(panelGeomX(end,k+1)-os*b))),...
                        -1*(panelGeomX(end,k+1)-os*b-VcontrolPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i),0);
                    
                    [u3,v3,w3] = VORTEX((-1*(panelGeomX(end,k)-os*b))- wakeT*b ,tan(beta)*((-1*(panelGeomX(end,k)-os*b))- wakeT*b),tan(alpha(d2))*((-1*(panelGeomX(end,k)-os*b))- wakeT*b),...
                        wakeT*b-(-1*VcontrolPointsX(j,i)),panelGeomY(k)-tan(beta)*((-1*(panelGeomX(end,k)-os*b))- wakeT*b)-controlPointsY(i),tan(alpha(d2))*(wakeT*b-(-1*panelGeomX(end,k))));

                    %Trailing Vorticies on the wing
                    %Vortex 4: quarter chord to TE
                    %Vortex 5: TE to quarter chord
                    [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)-os*b) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-VcontrolPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i),0);
                    [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k)-os*b)),0,0,(-1*(panelGeomX(end,k)-os*b))-(-1*VcontrolPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                    
                    %Sum contributions from each segment in the u,v,w
                    %directions
                    w = w1 + w2 + w3 + w4 + w5;
                    
                    %Apply effects of camber and dihedral to contributions
                    %before assigning to vorticity matrix. Note
                    %asymmetircal dihedral angle on each wing
                   if panelCount > N/2
                    wi = wi + gammaMatrix(m,k)*w*cos(phi)*cos(atan(camberMat(j,i)));
                   else
                    wi = wi + gammaMatrix(m,k)*w*cos(phi)*cos(atan(camberMat(j,i)));
                   end        
            end
        end 
       alphai = atan(-wi/Uinf);
       inducedDragi = potentialLiftDist(j,i)*sin(alphai);
       potIndDrag(j,i) = inducedDragi;
       if j <= NLE
           %Calculate Lift due to suction if panels part of LE
           normalForce(j,i) = potentialLiftDist(j,i)*cos(alpha(d2)) + inducedDragi*sin(alpha(d2));
           dxLE = panelGeomX(j,i+1) - panelGeomX(j,i);
           dy = panelGeomY(i+1) - panelGeomY(i);
           angleLE = atan(dxLE/dy);
           suctionNormalForce(j,i) = (normalForce(j,i)*sin(alpha(d2))- inducedDragi)/cos(angleLE);
       end
       panelCount = panelCount + 1;
    end
end
%Calculate Vortex Lift and Drag
vortexLift =(suctionNormalForce+normalForce)*cos(alpha(d2));
vortexInducedDrag = vortexLift*tan(alpha(d2));
%Calculate No Suction Lift and Drag
noSucLift = (normalForce)*cos(alpha(d2));
noSucInducedDrag = normalForce*tan(alpha(d2));
%Potential Lift and Drag on remaining planform
remainingPlanformLift = zeros(cPanels,2*sum(bPanels));
remainingPlanformLift(NLE+1:end,:) =  potentialLiftDist(NLE+1:end,:);
remainingPlanformDrag = zeros(cPanels,2*sum(bPanels));
remainingPlanformDrag(NLE+1:end,:) =  potIndDrag(NLE+1:end,:);
%Combine the remaining planform lift and drag with the leading edge Vortex
%Lift and Drag
combineLift = vortexLift + remainingPlanformLift;
combineDrag = vortexInducedDrag + remainingPlanformDrag;
%Combine the remaining planform lift and drag with the leading edge No
%Suction Lift and Drag
combineDragNV = noSucInducedDrag + remainingPlanformDrag;
combineLiftNV = noSucLift + remainingPlanformLift;


LiftPot = sum(sum(potentialLiftDist));
CLpot(d2) = 2*LiftPot/(rho*S*Uinf^2);
PotentialIndDrag = sum(sum(potIndDrag));
CDiPot(d2) = 2*PotentialIndDrag/(rho*S*Uinf^2);

TotalLift = sum(sum(combineLift));
CL(d2) = 2*TotalLift/(rho*S*Uinf^2);
TotalDrag = sum(sum(combineDrag));
CDi(d2) = 2*TotalDrag/(rho*S*Uinf^2);

end

%Function to calculate vorticity contributions

function [u, v, w] = VORTEX(Gx, Gy, Gz, Rx, Ry, Rz)
x1n = 0;
y1n = 0;
z1n = 0;

x2n = Gx;
y2n = Gy;
z2n = Gz;

%Set Control Point location

%Convention 1: Textbook convention(Control Point minus Vortex Root)
% x = Rx;
% y = Ry;
% z = Rz;

%Convention 2: Kroo convention(Vortex Root minus Control Point)
x = -Rx;
y = -Ry;
z = -Rz;



R1XR2x = 1*((y-y1n)*(z-z2n)-(y-y2n)*(z-z1n));
R1XR2y = -1*((x-x1n)*(z-z2n)-(x-x2n)*(z-z1n));
R1XR2z = 1*((x-x1n)*(y-y2n)-(x-x2n)*(y-y1n));

R1XR2 = R1XR2x^2 + R1XR2y^2 + R1XR2z^2;


sig1 = ( (x2n-x1n)*(x-x1n) + (y2n-y1n)*(y-y1n) + (z2n-z1n)*(z-z1n))/...
    sqrt((x-x1n)^2 + (y-y1n)^2 + (z-z1n)^2);
sig2 = ( (x2n-x1n)*(x-x2n) + (y2n-y1n)*(y-y2n) + (z2n-z1n)*(z-z2n))/...
    sqrt((x-x2n)^2 + (y-y2n)^2 + (z-z2n)^2);


TOL = 1.0E-10;
G2=Gx*Gx+Gy*Gy+Gz*Gz;
R2=Rx*Rx+Ry*Ry+Rz*Rz;
TOL2 = TOL * G2;
 
GXRx=Gy*Rz-Gz*Ry;
GXRy=Gz*Rx-Gx*Rz;
GXRz=Gx*Ry-Gy*Rx;

E1=GXRx*GXRx+GXRy*GXRy+GXRz*GXRz;
E2=(Gx+Rx)*(Gx+Rx)+(Gy+Ry)*(Gy+Ry)+(Gz+Rz)*(Gz+Rz);
ERROR = 0;
if R2 <= TOL2 ERROR = 1; end
if E2 <= TOL2 ERROR = 1; end
if E1 <= TOL2*R2 ERROR = 1; end

if ERROR==0
    V = (sig1-sig2) / (4.0*pi);
    u = V*(R1XR2x/R1XR2);
    v = V*(R1XR2y/R1XR2);
    w = V*(R1XR2z/R1XR2);
else
    u=0.0;
    v=0.0;
    w=0.0;
end
 


end

end

