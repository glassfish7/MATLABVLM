function [CLpot,CLNvortex,CLvortex,CDiPot,CDiNV,CDiTotal] = vlm3aSweep(sections,data,bPanels,cPanels,alpha,Uinf,rho,LE)
%VLM2SweepAlpha Sweeps Alpha for a given planform with VLM2
%   data = [taper;rootC;AR;LambLE] column for each section


%Set number of panels and alpha
NLE = LE; %Number of rows of panels to apply suction analogy to starting with the Leading Edge


%Use Geometry Engine to Obtain Panels
[controlPointsY,controlPointsX,panelQuarterC,panelGeomY,panelGeomX,VcontrolPointsX,S,b,N,K] = geomEngine(sections,data,bPanels,cPanels);

realAR = b^2/S;
fprintf('The Aspect Ratio of this Wing is :%f\n',realAR);
fprintf('The Wing Area is: %f\n',S);

CLvortex = zeros(1,length(alpha));
CLNvortex = zeros(1,length(alpha));
CLpot = zeros(1,length(alpha));
CDiPot = zeros(1,length(alpha));
CDiTotal = zeros(1,length(alpha));
CDiNV = zeros(1,length(alpha));


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
                
                   vorticityMatrix(panelCount,contCount) = ...
                   VXYZ(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),-1*(panelQuarterC(m,k)-controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i))...    
                    + VXYZ(400*b - (-1*panelQuarterC(m,k+1)),0,-1*(panelQuarterC(m,k+1)-controlPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i))...
                    + VXYZ((-1*panelQuarterC(m,k))- 400*b ,0,400*b-(-1*controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i))  ;
               
                   contCount = contCount +1;
                
                
            end
        end
        
        %Update counters
        panelCount = panelCount + 1;
        contCount =1;
    end
end

for z = 1:length(alpha)

%Generate RHS of matrix equation with no penetration boundary condition

RHS = -Uinf*sin(alpha(z))*ones(N,1);

%Solve for vorticity

%Using linsolve
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



sPoints = zeros(2,sections);
ePoints = zeros(2,sections);

for i = 1:sections
    
    sPoints(1,i) = 1 + sum(bPanels(sections:-1:i))-bPanels(i);
    sPoints(2,i) = sPoints(1,i) + sum(bPanels(1:i))+sum(bPanels(1:i-1));


    ePoints(1,i) = sPoints(1,i) + bPanels(i)-1;
    ePoints(2,i) = sPoints(2,i) + bPanels(i)-1;


end



for t = 1:sections
          potentialLiftDist(:,[sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])  = ...
              rho*Uinf*gammaMatrix(:,[sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])*K(t);

end



%Vortex Lift Section

%Calculate Zero Leading edge suction lift 
%Find Normal Force on each panel for zero-suction lift

normalForce = zeros(cPanels,2*sum(bPanels));
potIndDrag =  zeros(cPanels,2*sum(bPanels));
for i = 1:length(controlPointsY)
    for j = 1:cPanels
        wi = 0;
        %Calculate downwash on panel      
        for k = 1:length(controlPointsY)
            for m = 1:cPanels
                
                    cont = ...
                    VXYZ(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),-1*(panelQuarterC(m,k)-VcontrolPointsX(j,i)),panelGeomY(k)-controlPointsY(i))...    
                   +  VXYZ(400*b - (-1*panelQuarterC(m,k+1)),0,-1*(panelQuarterC(m,k+1)-VcontrolPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i))...
                     + VXYZ((-1*panelQuarterC(m,k))- 400*b ,0,400*b-(-1*VcontrolPointsX(j,i)),panelGeomY(k)-controlPointsY(i))  ;
                    
                 
                    wi = wi + gammaMatrix(m,k)*cont;
     
            end
        end
       
       alphai = atan(-wi/Uinf);
       gammaIncKrhoU = potentialLiftDist(j,i);
       inducedDragi = gammaIncKrhoU*sin(alphai);
       potIndDrag(j,i) = inducedDragi;
       
       if j <= NLE
        normalForce(j,i) = gammaIncKrhoU*cos(alpha(z)) + inducedDragi*sin(alpha(z));
       end
    end
end


%Calculate Lift due to suction
suctionNormalForce = zeros(cPanels,2*sum(bPanels));

for i = 1:length(controlPointsY)
    for j = 1:NLE
        wi = 0;
        
        %Calculate downwash on quarter chord position       
        for k = 1:length(controlPointsY)
            for m = 1:cPanels
                
                    cont = ...
                    VXYZ(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),-1*(panelQuarterC(m,k)-VcontrolPointsX(j,i)),panelGeomY(k)-controlPointsY(i))...    
                   +  VXYZ(400*b - (-1*panelQuarterC(m,k+1)),0,-1*(panelQuarterC(m,k+1)-VcontrolPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i))...
                     + VXYZ((-1*panelQuarterC(m,k))- 400*b ,0,400*b-(-1*VcontrolPointsX(j,i)),panelGeomY(k)-controlPointsY(i))  ;
                    
                 
                    wi = wi + gammaMatrix(m,k)*cont;
            end
        end
       
      
       dxLE = panelGeomX(j,i+1) - panelGeomX(j,i);
       dy = panelGeomY(i+1) - panelGeomY(i);
       angleLE = atan(dxLE/dy);
       gammaIncKrhoU = potentialLiftDist(j,i);
       alphai = atan(-wi/Uinf);
     
       inducedDragi = gammaIncKrhoU*sin(alphai);
       suctionNormalForce(j,i) = (normalForce(j,i)*sin(alpha(z))- inducedDragi)/cos(angleLE);
    end
end

vortexLift =(suctionNormalForce+normalForce)*cos(alpha(z));
noVortexLift = (normalForce)*cos(alpha(z));
vortexInducedDrag = vortexLift*tan(alpha(z));
noSucInducedDrag = normalForce*tan(alpha(z));


remainingPlanformLift = zeros(cPanels,2*sum(bPanels));
remainingPlanformLift(NLE+1:end,:) =  potentialLiftDist(NLE+1:end,:);

remainingPlanformDrag = zeros(cPanels,2*sum(bPanels));
remainingPlanformDrag(NLE+1:end,:) =  potIndDrag(NLE+1:end,:);


combineLift = vortexLift + remainingPlanformLift;

combineDrag = vortexInducedDrag + remainingPlanformDrag;
combineDragNV = noSucInducedDrag + remainingPlanformDrag;
combineLiftNV = noVortexLift + remainingPlanformLift;
%Kill Fuselage Lift if necessary

% if fuselage == 1
%     [~,cut] = size(controlPointsY);
%     combineLift(:,[cut/2 cut/2+1]) = 0;
%     potentialLiftDist(:,[cut/2 cut/2+1]) = 0;
% end


%Calculate spanwise circulation distribution by summing each column of the
%gamma matrix
totalLiftDistribution = sum(combineLift);


LiftPot = sum(sum(potentialLiftDist));
CLpot(z) = 2*LiftPot/(rho*S*Uinf^2);

% fprintf('CL Potential Lift = %f\n',CLPot);


% TotalLift = trapz([-b/2 controlPointsY b/2],[0 totalLiftDistribution 0]);
TotalLift = sum(totalLiftDistribution);
CLvortex(z) = 2*TotalLift/(rho*S*Uinf^2);
% fprintf('CL Vortex Lift = %f\n',CL);




NVLiftDistribution = sum(combineLiftNV);
NVTotalLift = sum(NVLiftDistribution);
CLNvortex(z) = 2*NVTotalLift/(rho*S*Uinf^2);
%fprintf('CL No Vortex Lift = %f\n',CLNV);





PotentialIndDrag = sum(sum(potIndDrag));

CDiPot(z) = 2*PotentialIndDrag/(rho*S*Uinf^2);
% fprintf('CDi = %f\n',CDiPot);

IndDrag = sum(sum(combineDrag));

CDiTotal(z) = 2*IndDrag/(rho*S*Uinf^2);
% fprintf('CDi Vortex = %f\n',CDiTotal);


NVIndDrag = sum(sum(combineDragNV));

CDiNV(z) = 2*NVIndDrag/(rho*S*Uinf^2);
%fprintf('CDi No Vortex = %f\n',CDiNV);





end

%Function to calculate vorticity contributions


function [VZ] = VXYZ (GX, GY, RX, RY)
    % Induced velocity in Z direction due to unit strength vortex:
    % GX is length of vortex in x direction;
    % RX is distance from vortex root to control point in the x direction;
    % GX=XT-XR; GY=YT-YR; RX=XR-XC; RY=YR-YC;
    TOL=1.0E-10; colinear=false; VZ=0.0;
    R2=RX*RX+RY*RY; G2=GX*GX+GY*GY;
    TOL2=TOL*G2; GXRZ=GX*RY-GY*RX;
    % Check to see if control point lies in line with vortex:
    E1=GXRZ*GXRZ; GR=GX*RX+GY*RY;
    E2=(GX+RX)^2+(GY+RY)^2;
    if(R2<=TOL2 || E2<=TOL2) colinear=true; end
    if(E1<=TOL2*R2) colinear=true; end
    if (~colinear) VZ=GXRZ*(GR/sqrt(R2)-(G2+GR)/sqrt(E2))/(4.0*pi*E1); end
end

function [Vx, Vy, Vz] = VortexVxyz(Gx, Gy, Gz, Rx, Ry, Rz)
    % Compute velocity induced by a unit strength vortex filament.
    % Vx, Vy, Vz are the velocity components in the coordinate directions.
    % Gx, Gy, Gz are the vortex lengths in each direction: tip - root.
    % Rx, Ry, Rz is the radius vector from vortex root to control point.
    

    TOL = 1.0E-10;
    R2=Rx*Rx+Ry*Ry+Rz*Rz;
    G2=Gx*Gx+Gy*Gy+Gz*Gz;
    GR=Gx*Rx+Gy*Ry+Gz*Rz;
    TOL2 = TOL * G2;
 

    GXRx=Gy*Rz-Gz*Ry;
    GXRy=Gz*Rx-Gx*Rz;
    GXRz=Gx*Ry-Gy*Rx;
 

    % Check to see if control point lies in line with vortex:
    E1=GXRx*GXRx+GXRy*GXRy+GXRz*GXRz;
    E2=(Gx+Rx)*(Gx+Rx)+(Gy+Ry)*(Gy+Ry)+(Gz+Rz)*(Gz+Rz);
    ERROR = 0;
    if R2 <= TOL2 ERROR = 1; end
    if E2 <= TOL2 ERROR = 1; end
    if E1 <= TOL2*R2 ERROR = 1; end
 

    if ERROR==0
       V = (GR/sqrt(R2) - (G2+GR)/sqrt(E2)) / (4.0*pi*E1);
       Vx = V*GXRx;
       Vy = V*GXRy;
       Vz = V*GXRz;
    else
       Vx=0.0;
       Vy=0.0;
       Vz=0.0;
    end
 

end

end

