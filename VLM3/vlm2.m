%2D VLM Vortex Lift
% Safa Bakhshi

%Set Planforms 
%   data = [taper;rootC;AR;LambLE(degrees)] column for each section
% Please un-comment the desired planform to run the code

%Simple Delta
% sections =1;
% taper = 0; %Taper Ratio
% S = 100; % Reference Area
% ARfix = 2; %Aspect Ratio
% b = sqrt(ARfix*S); %Span
% rootC = S*2/b; %Root Chord
% leLambda= (90-atand((b/2)/rootC)); %Sweep angle 
% bPanels = 16;
% cPanels =16;
% data = [taper;rootC;ARfix;leLambda];


%Arrow Wings
% sections =1;
% lambda = 70;
% AR1 = 1.03992;
% AR2 = 4*tand(90-lambda);
% AR3 = 1.99436;
% RC = 10;
% bPanels = 16;
% cPanels =16;

% -.4
% data = [0;RC-(-.4*RC);AR1;lambda];
% 0
% data = [0;RC;AR2;lambda];
% .27
% data = [0;RC-(.27*RC);AR3;lambda];


%Clipped Delta Wing
% sections =1;
% lambda = 63;
% AR = .873471;
% bPanels = 16;
% cPanels =16;
% taper = 0.4;
% data = [0.4;10;AR;lambda];


%Double Delta 75/65
% sections = 2;
% bPanels = [6,7];
% cPanels = 16;
% data = [0.39319,0;14.4533,5.6829;.466822,1.86524;75,65];



%765-072B
% sections =4;
% data = [0.648294,.698608,.345206,.205671;187.23,121.37,84.79,29.27;.075305,.089833,.573734,3.46289;84.9578,82.3695,68,40];
% bPanels = [3,2,8,15];
% cPanels = 16;
% fuselage = 0;



 

%Set number of panels and alpha
NLE = 2; %Number of rows of panels to apply suction analogy to starting with the Leading Edge
alpha = 25*(pi/180); % Angle of attack in degrees/converted to Radians
Uinf =1;
rho =1;

%Use Geometry Engine to Obtain Panels
[controlPointsY,controlPointsX,panelQuarterC,panelGeomY,panelGeomX,VcontrolPointsX,S,b,N,K] = geomEngine(sections,data,bPanels,cPanels);

realAR = b^2/S;
fprintf('The Aspect Ratio of this Wing is :%f\n',realAR);
fprintf('The Wing Area is: %f\n',S);
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


%Generate RHS of matrix equation with no penetration boundary condition

RHS = -Uinf*sin(alpha)*ones(N,1);

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

for i = 1:length(controlPointsY)
    for j = 1:NLE
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
       
       normalForce(j,i) = gammaIncKrhoU*cos(alpha) + inducedDragi*sin(alpha);
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
       suctionNormalForce(j,i) = (normalForce(j,i)*sin(alpha)- inducedDragi)/cos(angleLE);
    end
end

vortexLift =(suctionNormalForce+normalForce)*cos(alpha);

remainingPlanformLift = zeros(cPanels,2*sum(bPanels));
remainingPlanformLift(NLE+1:end,:) =  potentialLiftDist(NLE+1:end,:);

combineLift = vortexLift + remainingPlanformLift;

%Kill Fuselage Lift if necessary

% if fuselage == 1
%     [~,cut] = size(controlPointsY);
%     combineLift(:,[cut/2 cut/2+1]) = 0;
%     potentialLiftDist(:,[cut/2 cut/2+1]) = 0;
% end


%Plot Panel Geometry and Control Points
%Create a Y coordinate matrix for the contour plot
bigY = zeros(cPanels,2*sum(bPanels));
for n = 1:cPanels
    bigY(n,:) = controlPointsY;
end

%Plot a potential lift plot
figure(3)
contourf(bigY,VcontrolPointsX,potentialLiftDist,90,'linecolor','none')
colorbar
title('Potential lift contour plot of planform');
xlabel('y');
ylabel('-x');
colormap jet


%Plot a vortex lift contour plot
figure(4)
contourf(bigY,VcontrolPointsX,vortexLift,90,'linecolor','none')
colorbar
title('Vortex lift distribution of planform');
xlabel('y');
ylabel('-x');
colormap jet

%Plot a total lift contour plot
figure(5)
contourf(bigY,VcontrolPointsX,combineLift,90,'linecolor','none')
colorbar
title('Total lift distribution of planform');
xlabel('y');
ylabel('-x');
colormap jet


%Calculate spanwise Lift distribution by summing each column of the
%gamma matrix
potliftDistribution = sum(potentialLiftDist);
% Plot spanwise gamma distribution
figure(7)
plot([-b/2 controlPointsY b/2],[0 potliftDistribution 0])
title('Calculated spanwise Potential Flow Lift distribution');
xlabel('y')
ylabel('L''');

% LiftPot = rho*Uinf*trapz([-b/2 controlPointsY b/2],[0 potliftDistribution 0]);
LiftPot = sum(potliftDistribution);
CLPot = 2*LiftPot/(rho*S*Uinf^2);

fprintf('CL Potential Lift = %f\n',CLPot);


%Calculate spanwise circulation distribution by summing each column of the
%gamma matrix
totalLiftDistribution = sum(combineLift);

% Plot spanwise gamma distribution
figure(8)
plot([-b/2 controlPointsY b/2],[0 totalLiftDistribution 0])
title('Calculated spanwise Vortex Lift distribution');
xlabel('y')
ylabel('L''');

% TotalLift = trapz([-b/2 controlPointsY b/2],[0 totalLiftDistribution 0]);
TotalLift = sum(totalLiftDistribution);

CL = 2*TotalLift/(rho*S*Uinf^2);

fprintf('CL Vortex Lift = %f\n',CL);








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
