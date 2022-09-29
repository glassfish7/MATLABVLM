

AR = linspace(0.5,3,6);
CLAR = [0.426963, 0.72, 1.006365, 1.245282, 1.438087,1.621500];

alpha  = linspace(0,30,7);
CLALPHA = [0,0.177469, 0.343242, 0.541046, 0.741104, 0.936342, 1.133304];


figure(1)
hold on
title('Aspect Ratio(AR) vs C_{L}','interpreter','tex')
xlabel('AR')
ylabel('C_L')
plot(AR,CLAR)
scatter(AR,CLAR)
grid on

figure(2)
hold on
title('\alpha vs C_{L}','interpreter','tex')
xlabel('\alpha')
ylabel('C_L')
plot(alpha,CLALPHA)
scatter(alpha,CLALPHA)
grid on
