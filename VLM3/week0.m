% AA Research Week Zero Driver Script

%Simple Delta
sections =1;
taper = 0; %Taper Ratio
S = 100; % Reference Area
ARfix = .25; %Aspect Ratio
b = sqrt(ARfix*S); %Span
rootC = S*2/b; %Root Chord
leLambda= (90-atand((b/2)/rootC)); %Sweep angle 
bPanels = 16;
cPanels =16;
data = [taper;rootC;ARfix;leLambda];
alpha = linspace(0,14,15)*(pi/180);
Uinf =1;
rho =1;
LE =4;

%No longer works: see later code
[CLvortex,~,~,CDiTotal] = vlm3aSweep(sections,data,bPanels,cPanels,alpha,Uinf,rho,LE);

alpha = linspace(0,36,35)*(pi/180);

[~,CLpot,CDiPot,~] = vlm3aSweep(sections,data,bPanels,cPanels,alpha,Uinf,rho,LE);

figure(3)
hold on
plot(CLpot.^2,CDiPot);
plot(CLvortex.^2,CDiTotal);
title('Planform Drag Polar')
xlabel('C_L^2')
ylabel('C_{Di}');
grid on
legend('Potential','Vortex','Location','best');