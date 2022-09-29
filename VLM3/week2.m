% AA Research Week Two Driver Script

%Simple Delta
sections =1;
taper = 0; %Taper Ratio
S = 100; % Reference Area
ARfix = 2; %Aspect Ratio
b = sqrt(ARfix*S); %Span
rootC = S*2/b; %Root Chord
leLambda= (90-atand((b/2)/rootC)); %Sweep angle 
bPanels = 16;
cPanels =16;
data = [taper;rootC;ARfix;leLambda];
alpha = linspace(0,14,15)*(pi/180);
Uinf =1;
rho =1;
LE =2;


[~,~,CLvortex,~,~,CDiTotal] = vlm3aSweep(sections,data,bPanels,cPanels,alpha,Uinf,rho,LE);

alpha = linspace(0,18,19)*(pi/180);

[CLpot,~,~,CDiPot,~,~] = vlm3aSweep(sections,data,bPanels,cPanels,alpha,Uinf,rho,LE);

alpha = linspace(0,20,21)*(pi/180);

[~,CLNvortex,~,~,CDiNV,~] = vlm3aSweep(sections,data,bPanels,cPanels,alpha,Uinf,rho,LE);


figure(4)
hold on
plot(CLpot.^2,CDiPot);
plot(CLvortex.^2,CDiTotal);
plot(CLNvortex.^2,CDiNV);
title('Planform Drag Polar')
xlabel('C_L^2')
ylabel('C_{Di}');
grid on
legend('Potential','Vortex','No Vortex','Location','best');