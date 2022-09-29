%Safa Bakhshi

%Week 10 Presentation


%% 765-072B vehicle

sections =4;
Uinf =1;
rho =1;
LE =2;
alpha  = linspace(0,25,26);
data = [0.648294,.698608,.345206,.205671;187.23,121.37,84.79,29.27;.075305,.089833,.573734,3.46289;84.9578,82.3695,68,40];
bPanels = [3,2,8,15];
cPanels = 16;
fuselage =0;

[CLvortex,CLpot] = vlm2aSweep(sections,data,bPanels,cPanels,alpha*(pi/180),Uinf,rho,LE,fuselage);


figure()
hold on
plot(alpha,CLvortex);
plot(alpha,CLpot);
title('C_L vs \alpha for 765-072B Vehicle','interpreter','tex');
xlabel('\alpha(degrees)');
ylabel('C_L');
legend('Vortex Lift','Potential','Location','best')
grid on

