%Safa Bakhshi
%Code for Week 5

%% Slender Delta Wing vs Arrow Wing

%Slender Delta Wing Alpha 

%Define Planform
rootC = 30; %Root Chord
taper = 0; %Taper Ratio
S = 100; % Reference Area
ARfix = 0.5; %Aspect Ratio
b = sqrt(ARfix*S); %Span
leLambda= (90-atand((b/2)/rootC))*(pi/180); %Sweep angle in degrees/converted to Radians

%Set number of panels and alpha
N = 16;  %Number of panels = 2*N^2
alpha = linspace(0,20,21)*(pi/180); % Angle of attack in degrees/converted to Radians
Uinf =1;
rho =1;

CLVLM1 = vlmAlphaSweep(rootC,taper,S,ARfix,leLambda,N,alpha,Uinf,rho);

%Arrow Wing
leLambda = 85*(pi/180);

CLVLM2 = vlmAlphaSweep(rootC,taper,S,ARfix,leLambda,N,alpha,Uinf,rho);


%Slender Delta Wing Alpha

alphaFix =5*(pi/180);
AR = [.25:.25:2];
delta = 1;
leLambda= (90-atand((b/2)/rootC))*(pi/180); %Sweep angle in degrees/converted to Radians

CLVLM3 = vlmARSweep(rootC,taper,S,AR,leLambda,N,alphaFix,Uinf,rho,delta);

%Arrow Wing
leLambda = 85*(pi/180);
delta = 0;
CLVLM4 = vlmARSweep(rootC,taper,S,AR,leLambda,N,alphaFix,Uinf,rho,delta);



CLslender1 = (ARfix*pi*sin(alpha))/2;

CLslender2 = (AR*pi*sin(alphaFix))/2;


figure(2)
hold on
plot(alpha,CLVLM1)
plot(alpha,CLVLM2)
plot(alpha,CLslender1)
legend('Delta','Arrow','Slender','Location','best');
title('C_{L} vs \alpha for Delta and Arrow Wings','interpreter','tex')
xlabel('\alpha')
ylabel('C_{L}');
grid on


figure(3)
hold on
plot(AR,CLVLM3)
plot(AR,CLVLM4)
plot(AR,CLslender2)
legend('Delta','Arrow','Slender','Location','best');
title('C_{L} vs AR for  for Delta and Arrow Wings','interpreter','tex')
xlabel('AR')
ylabel('C_{L}');
grid on