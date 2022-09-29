%Safa Bakhshi

%Week 9 Presentation

%% Arrow Wing
LE = 2;

alpha  = linspace(0,25,26);

%3 Notch ratios

sections =1;
lambda = 70;
AR1 = 1.03992;
AR2 = 4*tand(90-lambda);
AR3 = 1.99436;
RC = 10;
bPanels = 16;
cPanels =16;
Uinf =1;
rho =1;
% -.4

data1 = [0;RC-(-.4*RC);AR1;lambda];

% 0

data2 = [0;RC;AR2;lambda];

% .27


data3 = [0;RC-(.27*RC);AR3;lambda];


[CLvortex1,CLpot1] = vlm2aSweep(sections,data1,bPanels,cPanels,alpha*(pi/180),Uinf,rho,LE);

[CLvortex2,CLpot2] = vlm2aSweep(sections,data2,bPanels,cPanels,alpha*(pi/180),Uinf,rho,LE);

[CLvortex3,CLpot3] = vlm2aSweep(sections,data3,bPanels,cPanels,alpha*(pi/180),Uinf,rho,LE);

figure()
hold on
plot(alpha,CLvortex1);
plot(alpha,CLpot1);
title('C_{L} vs \alpha for Arrow Wing a/X = -.4');
xlabel('\alpha(degrees)');
ylabel('C_{L}');
legend('Vortex Lift','Potential','Location','best')
grid on

figure()
hold on
plot(alpha,CLvortex2);
plot(alpha,CLpot2);
title('C_{L} vs \alpha  for Arrow Wing a/X = 0');
xlabel('\alpha(degrees)');
ylabel('C_{L}');
legend('Vortex Lift','Potential','Location','best')
grid on

figure()
hold on
plot(alpha,CLvortex3);
plot(alpha,CLpot3);
title('C_{L} vs \alpha for Arrow Wing a/X = .27','Interpreter','tex');
xlabel('\alpha(degrees)');
ylabel('C_{L}');
legend('Vortex Lift','Potential','Location','best')
grid on

figure()
hold on
plot([-.4 0 .27],[CLvortex1(end) CLvortex2(end) CLvortex3(end)]);
title('C_{L} vs a/X for Arrow Wings with \Lambda = 70','Interpreter','tex');
xlabel('a/X');
ylabel('C_{L}');
legend('Vortex Lift','Location','best')
grid on


%% Clipped Delta

sections =1;
lambda = 63;
AR = .873471;
bPanels = 16;
cPanels =16;
Uinf =1;
rho =1;
taper = 0.4;
LE =2;
alpha  = linspace(0,35,36);


data4 = [0.4;10;AR;lambda];

[CLvortex4,CLpot4] = vlm2aSweep(sections,data4,bPanels,cPanels,alpha*(pi/180),Uinf,rho,LE);

figure()
hold on
plot(alpha,CLvortex4);
plot(alpha,CLpot4);
title('C_{L} vs \alpha for Clipped Delta \lambda=0.4','Interpreter','tex');
xlabel('\alpha(degrees)');
ylabel('C_{L}');
legend('Vortex Lift','Potential','Location','best')
grid on


%% Double Delta Wing

sections = 2;
Uinf =1;
rho =1;
LE =2;
alpha  = linspace(0,25,26);
bPanels = [6,7];
cPanels = 8;
data5 = [0.39319,0;14.4533,5.6829;.466822,1.86524;75,65];


[CLvortex5,CLpot5] = vlm2aSweep(sections,data5,bPanels,cPanels,alpha*(pi/180),Uinf,rho,LE);

figure()
hold on
plot(alpha,CLvortex5);
plot(alpha,CLpot5);
title('C_{L} vs \alpha for 75/65 degree Double Delta Wing','Interpreter','tex');
xlabel('\alpha(degrees)');
ylabel('C_L');
legend('Vortex Lift','Potential','Location','best')
grid on
