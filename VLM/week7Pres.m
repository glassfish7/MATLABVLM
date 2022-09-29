%Safa Bakhshi
%Code for Week 7

%% Slender Delta Wing

%Slender Delta Wing Alpha 

%Define Planform
taper = 0; %Taper Ratio
S = 100; % Reference Area
ARfix = 2; %Aspect Ratio
b = sqrt(ARfix*S); %Span
rootC = S*2/b; %Root Chord
leLambda= (90-atand((b/2)/rootC))*(pi/180); %Sweep angle in degrees/converted to Radians

%Set number of panels and alpha
N = 16;  %Number of panels = 2*N^2
alpha = linspace(0,25,26)*(pi/180); % Angle of attack in degrees/converted to Radians
Uinf =1;
rho =1;

CLVLM1 = vlmAlphaSweep(rootC,taper,S,ARfix,leLambda,N,alpha,Uinf,rho);
CLVLM2 = vlmAlphaSweepSuc(rootC,taper,S,ARfix,leLambda,N,alpha,Uinf,rho);

%Slender Delta Wing Alpha

alphaFix =15*(pi/180);
AR = [.25:.25:2];
delta = 1;

CLVLM3 = vlmARSweep(rootC,taper,S,AR,leLambda,N,alphaFix,Uinf,rho,delta);
CLVLM4 = vlmARSweepSuc(rootC,taper,S,AR,leLambda,N,alphaFix,Uinf,rho,delta);

CLslender1 = (ARfix*pi*sin(alpha))/2;

CLslender2 = (AR*pi*sin(alphaFix))/2;

%NLE Sweep
NLE = [1:N];
alphaFix =25*(pi/180);
CLVLM5 = vlmNLESweepSuc(rootC,taper,S,ARfix,leLambda,N,alphaFix,Uinf,rho,delta,NLE);

ARfix = 0.5;

CLVLM6 = vlmNLESweepSuc(rootC,taper,S,ARfix,leLambda,N,alphaFix,Uinf,rho,delta,NLE);





figure(2)
hold on
plot(alpha*(180/pi),CLVLM1)
plot(alpha*(180/pi),CLVLM2)
plot(alpha*(180/pi),CLslender1)
legend('Potential VLM','Vortex Lift VLM','Slender Wing Theory','Location','best');
title('C_{L} vs \alpha for Slender Delta Wing','interpreter','tex')
xlabel('\alpha(degrees)')
ylabel('C_{L}');
grid on


figure(3)
hold on
plot(AR,CLVLM3)
plot(AR,CLVLM4)
plot(AR,CLslender2)
legend('Potential VLM','Vortex Lift VLM','Slender Wing Theory','Location','best');
title('C_{L} vs AR for Slender Delta Wing','interpreter','tex')
xlabel('AR')
ylabel('C_{L}');
grid on

figure(4)
hold on
plot(NLE,CLVLM5)
legend('VLM','Polhamus','Location','best');
title('C_{L} vs NLE for AR =2','interpreter','tex')
xlabel('NLE')
ylabel('C_{L}');
grid on


figure(5)
hold on
plot(NLE,CLVLM6)
legend('VLM','Polhamus','Location','best');
title('C_{L} vs NLE for AR =0.5','interpreter','tex')
xlabel('NLE')
ylabel('C_{L}');
grid on

