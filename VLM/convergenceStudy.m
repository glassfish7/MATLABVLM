%Dirty Plots

N1 = [32,128,512];
Clp1 = [0.949559,0.934721,0.930579];

figure(1)
plot(N1,Clp1);
ylim([0.5 1])
ylabel('C_{L_p}');
xlabel('NPanels');
title('Potential lift convergence study for simple delta wing');
grid on

N2 = [56,208, 823];
Cpl2 = [0.769366, 0.764854, 0.764654];

figure(2)
plot(N2,Cpl2);
ylim([0.5 1])
ylabel('C_{L_p}');
xlabel('NPanels');
title('Potential lift convergence study for 75/65 degrees double delta wing');
grid on

N3 = [896,1792];
Cpl3 = [0.948909, 0.949014];

figure(3)
plot(N3,Cpl3);
ylim([0.5 1])
ylabel('C_{L_p}');
xlabel('NPanels');
title('Potential lift convergence study convergence study for 765-072B');
grid on