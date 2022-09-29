sections = 3;
data = [.6,0.3333,.521327;1067,753,248.49;.117151,0.845766,2.35319;80.95,66,42];
% sections = 1;
% data = [.6;1067;.083126;85.2];
bPanels = [3,8,8];
cPanels = 16;
camberLine = 0;
alphaIncRoot = 0*(pi/180); %Incidence at Root converted to radians
alphaIncTip = 0*(pi/180); %Incidence at Tip converted to radians


LE = 2;

alpha = linspace(0,35,36)*(pi/180);
Uinf = 1;
rho = 1;


[~,~,CLvortex,~,~,CDiTotal]= vlm3aSweep(sections,data,bPanels,cPanels,alpha,Uinf,rho,LE);


figure()
hold on
plot(alpha,CLvortex);
title('Aircraft Lift Curve')
xlabel('\alpha')
ylabel('C_{L}');
grid on

figure()
hold on
plot(CDiTotal,CLvortex);
title('Planform Drag Polar')
xlabel('C_{Di}')
ylabel('C_L');
grid on