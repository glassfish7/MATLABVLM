%Safa Bakhshi

%Simple Rectangular
sections = 1;
taper = 1;
S = 100;
ARfix = 10;
b = sqrt(ARfix*S); %Span
rootC = S/b; %Root Chord
phi = 0*(pi/180);%Dihedral Angle
alphaIncRoot = 0*(pi/180); %Incidence at Root converted to radians
alphaIncTip = 0*(pi/180); %Incidence at Tip converted to radians
leLambda = 0;
data = [taper;rootC;ARfix;leLambda];
camberLine = 0;
% camberLine = '2412';
% camberLine = [0.0468750000000000,0.109375000000000,0.171875000000000,0.234375000000000,0.296875000000000,0.359375000000000,0.421875000000000,...
% 0.484375000000000,0.546875000000000,0.609375000000000,0.671875000000000,0.734375000000000,0.796875000000000,0.859375000000000,0.921875000000000,...
% 0.984375000000000;0.0130723741319444,0.0153076171875000,0.0171088324652778,0.0184760199652778,0.0194091796875000,0.0199083116319444,0.0199734157986111,...
% 0.0196044921875000,0.0188015407986111,0.0175645616319444,0.0158935546875000,0.0137885199652778,0.0112494574652778,0.00827636718750000,0.00486924913194445,...
% 0.00102810329861111];

%Simple Delta
% sections =1;
% taper = 0; %Taper Ratio
% S = 100; % Reference Area
% ARfix = 2; %Aspect Ratio
% b = sqrt(ARfix*S); %Span
% rootC = S*2/b; %Root Chord
% alphaIncRoot = 0*(pi/180); %Incidence at Root converted to radians
% alphaIncTip = 0*(pi/180); %Incidence at Tip converted to radians
% leLambda= (90-atand((b/2)/rootC)); %Sweep angle 
% phi = 0*(pi/180); %Dihedral Angle
% data = [taper;rootC;ARfix;leLambda];
% camberLine = 0;

%VLM Control: Set panels and modelling options
bPanels = 16; %Spanwise Panels/2
cPanels = 16; %Chordwise Panels
wakeT = 400; %span lengths for wake termination
os = 0; %Kink Trailing Edge Offset Factor

%Freestream Control: Set alpha,beta, and leading edge panels for suction
%analogy, Uinf, and rho
NLE = 0; %Number of rows of panels to apply suction analogy to starting with the Leading Edge
alpha = (0:15)*(pi/180); % Angle of attack in degrees/converted to Radians 
beta = 0*(pi/180); %Sideslip angle in degrees/converted to Radians
Uinf =1; %Freestream velocity 
rho =1; %Freestream Density


[CLpot,CDiPot,CL,CDi,CDiTreffz] = vlm5aSweep(sections,data,alphaIncRoot,alphaIncTip,phi,camberLine,bPanels,cPanels,alpha,beta,Uinf,rho,NLE,wakeT,os);


figure()
plot(alpha*(180/pi),CLpot);
title('\alpha vs CL Potential','interpreter','tex')
xlabel('\alpha')
ylabel('C_L')
grid on

figure()
plot(alpha*(180/pi),CDiPot);
title('\alpha vs CDi Potential','interpreter','tex')
xlabel('\alpha')
ylabel('C_D_{i}')
grid on

figure()
plot(alpha*(180/pi),CDiTreffz);
title('\alpha vs CDi Potential (Trefftz)','interpreter','tex')
xlabel('\alpha')
ylabel('C_D_{i}')
grid on