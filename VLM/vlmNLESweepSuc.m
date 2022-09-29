function [CL] = vlmNLESweepSuc(rootC,taper,S,AR,leLambda,N,alpha,Uinf,rho,delta,LE)


for r = 1:length(LE)

tipC = rootC*taper; % Tip Chord
b = sqrt(AR*S); %Span
NLE = LE(r); %Number of rows of panels to apply suction analogy to starting with the Leading Edge
if delta == 1
leLambda= (90-atand((b/2)/rootC))*(pi/180); %Enable line if Delta
rootC = S*2/b;
end

%Calculate Root and Tip ending and starting points 
rootEnd = 0;
rootStart = rootC + rootEnd;

tipStart = rootStart - (b/2)*tan(leLambda);
tipEnd = tipStart - tipC;

% %Plot Planform
% figure(1)
% plot([0 0 -b/2 -b/2 0 b/2 b/2 0],[rootStart rootEnd tipEnd tipStart rootStart tipStart tipEnd rootEnd]); 
% title('Planform Shape and Dimensions')
% xlabel('y');
% ylabel('x');


%Generate Panels

%Descretize wing root and wing tip into N+1 points  
rootChordPoints = rootStart:-(rootStart-rootEnd)/N:rootEnd;
tipChordPoints =  tipStart:-(tipStart-tipEnd)/N:tipEnd;

if isempty(tipChordPoints)
    tipChordPoints = tipStart*ones(1,length(rootChordPoints));
end


%Divide wing spanwise into 2N+1 stations
K = b/(2*N);
panelGeomY = -b/2:K:b/2;

%Find wing centerline index
centreIndex = find(panelGeomY==0);


%Create panel x coordiantes
panelGeomX  = zeros(length(rootChordPoints),length(panelGeomY));
for i = 1:length(rootChordPoints)
    
    %Left Wing
    panelGeomX(i,1:centreIndex) = rootChordPoints(i) + ((tipChordPoints(i)-rootChordPoints(i))/(-b/2))*(panelGeomY(1:centreIndex));
    %Right Wing
    panelGeomX(i,centreIndex+1:end) = rootChordPoints(i) + ((tipChordPoints(i)-rootChordPoints(i))/(b/2))*(panelGeomY(centreIndex+1:end));

end


%Create quarter chord points for each panel
panelQuarterC = zeros(N,length(panelGeomY));
for i = 1:length(panelGeomY)
    for j = 1:N
        panelQuarterC(j,i) = panelGeomX(j,i) + (panelGeomX(j+1,i)-panelGeomX(j,i))/4;
    end
end


%Create three quarter chord points for each panel
tquarterPointsX = zeros(N,length(panelGeomY));
for i = 1:length(panelGeomY)
    for j = 1:N
        tquarterPointsX(j,i) = panelGeomX(j,i) + 3*(panelGeomX(j+1,i)-panelGeomX(j,i))/4;
    end
end


%Create control point x coordinate for each panel
controlPointsX = zeros(N,length(panelGeomY)-1);
for i= 1:length(panelGeomY)-1
    for j = 1:N
        controlPointsX(j,i) = (tquarterPointsX(j,i+1)+tquarterPointsX(j,i))/2;
    end
end


%Create on vortex control points for vortex lift calculations
VcontrolPointsX = zeros(N,length(panelGeomY)-1);
for i= 1:length(panelGeomY)-1
    for j = 1:N
        VcontrolPointsX(j,i) = (panelQuarterC(j,i+1)+panelQuarterC(j,i))/2;
    end
end


%Create control point y coordinate for each panel
controlPointsY= zeros(1,length(panelGeomY)-1);

for i = 1:length(panelGeomY)-1
    controlPointsY(i) = panelGeomY(i) + (panelGeomY(i+1)-panelGeomY(i))/2;
end


%Generate Vorticity Matrix

%Panels are numbered sequentially from 1 to the total number of panels
%Each panel is numbered in Figure 2

%The row in the matrix corresponding to a particular panel contains all of
%the downwash contribution on that panel from every other panel in the
%system

%The column in the matrix corresponding to a particular panel contains all
%of the contribution of that particular panel to every other panel in the
%system


vorticityMatrix = zeros(2*N^2,2*N^2);

%Set counters to 1
panelCount =1;
contCount =1;


%Loop through each panel on the wing. Calculate downwash contributions ONTO this
%panel 
for i = 1:length(controlPointsY)
    for j = 1:N
        
    
        %Loop through each panel on the wing. Calculate downwash contributions BY this
        %panel 
        for k = 1:length(controlPointsY)
            for m = 1:N
                  
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

RHS = -Uinf*sin(alpha)*ones(2*N^2,1);

%Solve for vorticity

%Using linsolve
gamma = linsolve(vorticityMatrix,RHS);


%Place gamma vector back into array for easy interpretation of results

gammaMatrix = zeros(N,2*N);

gammaCount =1;
for p = 1:2*N
    for o = 1:N
        gammaMatrix(o,p) = gamma(gammaCount);
        gammaCount = gammaCount + 1;
    end
end

potentialLiftDist  = rho*Uinf*gammaMatrix*K;

%Vortex Lift Section

%Calculate Zero Leading edge suction lift 
%Find Normal Force on each panel for zero-suction lift

normalForce = zeros(N,2*N);

for i = 1:length(controlPointsY)
    for j = 1:NLE
        wi = 0;
        %Calculate downwash on panel      
        for k = 1:length(controlPointsY)
            for m = 1:N
                
                    cont = ...
                    VXYZ(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),-1*(panelQuarterC(m,k)-VcontrolPointsX(j,i)),panelGeomY(k)-controlPointsY(i))...    
                   +  VXYZ(400*b - (-1*panelQuarterC(m,k+1)),0,-1*(panelQuarterC(m,k+1)-VcontrolPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i))...
                     + VXYZ((-1*panelQuarterC(m,k))- 400*b ,0,400*b-(-1*VcontrolPointsX(j,i)),panelGeomY(k)-controlPointsY(i))  ;
                    
                 
                    wi = wi + gammaMatrix(m,k)*cont;
     
            end
        end
       
       alphai = atan(-wi/Uinf);
       gamma = gammaMatrix(j,i);
       inducedDragi = rho*Uinf*gamma*K*sin(alphai);
       
       normalForce(j,i) = rho*Uinf*gamma*K*cos(alpha) + inducedDragi*sin(alpha);
    end
end


%Calculate Lift due to suction
suctionNormalForce = zeros(N,2*N);

for i = 1:length(controlPointsY)
    for j = 1:NLE
        wi = 0;
        
        %Calculate downwash on quarter chord position       
        for k = 1:length(controlPointsY)
            for m = 1:N
                
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
       gamma = gammaMatrix(j,i);
       alphai = atan(-wi/Uinf);
     
       inducedDragi = rho*Uinf*gamma*K*sin(alphai);
       suctionNormalForce(j,i) = (normalForce(j,i)*sin(alpha)- inducedDragi)/cos(angleLE);
    end
end

vortexLift =(suctionNormalForce+normalForce)*cos(alpha);

remainingPlanformLift = zeros(N,2*N);
remainingPlanformLift(NLE+1:end,:) =  potentialLiftDist(NLE+1:end,:);

combineLift = vortexLift + remainingPlanformLift;


%Calculate spanwise circulation distribution by summing each column of the
%gamma matrix
liftDistribution = sum(combineLift);


% Lift = trapz([-b/2 controlPointsY b/2],[0 liftDistribution 0]);
Lift = sum(liftDistribution);

CL(r) = 2*Lift/(rho*S*Uinf^2);

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

