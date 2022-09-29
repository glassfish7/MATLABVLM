%Simple Conical
clear all

%Define Simple Delta Wing
sections =1;
taper = 0; %Taper Ratio
S = 100; % Reference Area
ARfix =1; %Aspect Ratio
b = sqrt(ARfix*S); %Span
rootC = S*2/b; %Root Chord
leLambda= (90-atand((b/2)/rootC)); %Sweep angle 
data = [taper;rootC;ARfix;leLambda];

%Freestream Control: Set alpha, Uinf, and rho
alpha = 20*(pi/180); % Angle of attack in degrees/converted to Radians 
Uinf =1; rho = 1; %Freestream velocity and Density 

%VLM Control: Set panels and modelling options
bPanels = 12; %Spanwise panels divided by 2
cPanels = 24; %Chordwise Panels
wakeT = 400; %span lengths for wake termination

%Use Geometry Engine to Obtain Panels and Control Point Geometry data
[controlPointsY,controlPointsX,panelQuarterC,panelQuarterCX,panelHalfCY,panelTQuarterC,panelGeomY,panelGeomX,VcontrolPointsX,VcontrolPointsY,S,b,N,K,CAve,geomChord,chord]...
    = geomEngine3Conc(sections,data,bPanels,cPanels,"NoPlot");

xw = zeros(1,2*bPanels+1); %Trailing edge vortex release points x coordinates
yw = panelGeomY(end,:);    %Trailing edge vortex release points y coordinates
zw = zeros(1,2*bPanels+1); %Trailing edge vortex release points z coordinates

%VLM Main Process
vorticityMatrix = zeros(N,N); %Initialize AIC Matrix
%Set counters to 1
panelCount =1; contCount =1;
%Loop Through all panels
for i = 1:2*bPanels
    for j = 1:cPanels
        for k = 1:2*bPanels
            for m = 1:cPanels
                %Contribution from quarter chord bound vortex
                [u1,v1,w1] = VORTEX(0,panelQuarterC(m,k+1)-panelQuarterC(m,k),0,-1*(panelQuarterCX(m,k)-controlPointsX(j,i)),panelQuarterC(m,k)-controlPointsY(j,i),0);
                %Contribution from trailing semi infinite vortices
                [u2inf,v2inf,w2inf] = VORTEX(wakeT*b - (-1*xw(k+1)),0,0,...
                    -1*(xw(k+1)-controlPointsX(j,i)),yw(k+1)-controlPointsY(j,i),0);
                [u3inf,v3inf,w3inf] = VORTEX((-1*xw(k))- wakeT*b ,0,0,...
                    wakeT*b-(-1*controlPointsX(j,i)),yw(k)-controlPointsY(j,i),0);
                %Contrubtion from trailing vortex segments on wings
                [u4,v4,w4] = VORTEX((-1*(panelGeomX(end)) - (-1*panelQuarterCX(m,k))),panelGeomY(end,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-controlPointsX(j,i)),panelQuarterC(m,k+1)-controlPointsY(j,i),0);
                [u5,v5,w5] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(end))),panelQuarterC(m,k)-panelGeomY(end,k),0,-1*(panelGeomX(end)-controlPointsX(j,i)),panelGeomY(end,k)-controlPointsY(j,i),0);
                w = w1 + w2inf + w3inf + w4 + w5; %Sum up downwash contributions from each segment

                vorticityMatrix(panelCount,contCount) =  w;
                contCount = contCount + 1; 
            end
        end
        panelCount = panelCount + 1; contCount = 1;
    end
end

RHS = (Uinf*(-sin(alpha))).*ones(N,1); %Set RHS
gamma = linsolve(vorticityMatrix,RHS); %Solve
gammaMatrix = reshape(gamma,cPanels,2*sum(bPanels)); %Arrange Gamma into a matrix

%Use Kutta-Joukowski Theorem to calculate lift distribution
Kwid = repmat(diff(panelGeomY(end,:)),cPanels,1);%Calculate panel widths at quarter chord
potentialLiftDist  = rho*Uinf*gammaMatrix.*Kwid;
LiftPot = sum(sum(potentialLiftDist)); 
CLPot = 2*LiftPot/(rho*S*Uinf^2);
fprintf('CL = %f\n',CLPot);


%Plot contour plot of the vorticity distribution    
figure()
contourf(VcontrolPointsY,VcontrolPointsX,gammaMatrix,90,'linecolor','none')
colorbar
title('Total gamma distribution of planform');
xlabel('y');
ylabel('-x');
colormap jet   


% Function to calculate 
function [u,v,w] = VORTEX(Gx, Gy, Gz, Rx, Ry, Rz)
%Vortex Root at Origin
x1n = 0;
y1n = 0;
z1n = 0;

%Vortex Tip Location
x2n = Gx;
y2n = Gy;
z2n = Gz;

% Set Control Point location

%Convention 2: Kroo convention(Vortex Root minus Control Point)
x = -Rx;
y = -Ry;
z = -Rz;

% 
G2=Gx*Gx+Gy*Gy+Gz*Gz;
R2=Rx*Rx+Ry*Ry+Rz*Rz;
E2=(Gx+Rx)*(Gx+Rx)+(Gy+Ry)*(Gy+Ry)+(Gz+Rz)*(Gz+Rz);
R1Length = sqrt(R2);
GLength = sqrt(G2);
R2Length = sqrt(E2);
spp = (R1Length+GLength+R2Length)/2;
h = 2*sqrt(spp*(spp-R1Length)*(spp-GLength)*(spp-R2Length))/GLength;

% h = norm(cross([Rx,Ry,Rz],[Gx,Gy,Gz]))/norm([Gx,Gy,Gz]);
RC = 20;
if h < (1E-6*RC)
% if h < sqrt(4*.001*timeGlobal)
    u=0.0;
    v=0.0;
    w=0.0;
else
    %R1 X R2
    R1XR2x = 1*((y-y1n)*(z-z2n)-(y-y2n)*(z-z1n));
    R1XR2y = -1*((x-x1n)*(z-z2n)-(x-x2n)*(z-z1n));
    R1XR2z = 1*((x-x1n)*(y-y2n)-(x-x2n)*(y-y1n));

    %|R1 X R2|
    R1XR2 = R1XR2x^2 + R1XR2y^2 + R1XR2z^2;

    %sigma: see Virginia Tech textbook ch 6 pg. 18
    sig1 = ( (x2n-x1n)*(x-x1n) + (y2n-y1n)*(y-y1n) + (z2n-z1n)*(z-z1n))/...
        sqrt((x-x1n)^2 + (y-y1n)^2 + (z-z1n)^2);
    sig2 = ( (x2n-x1n)*(x-x2n) + (y2n-y1n)*(y-y2n) + (z2n-z1n)*(z-z2n))/...
        sqrt((x-x2n)^2 + (y-y2n)^2 + (z-z2n)^2);

    V = (sig1-sig2) / (4.0*pi);
    u = V*(R1XR2x/R1XR2);
    v = V*(R1XR2y/R1XR2);
    w = V*(R1XR2z/R1XR2);
end

end
