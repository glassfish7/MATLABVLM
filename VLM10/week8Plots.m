%Quick Plots

AR = [0.5,1,1.5,2];
CL = [0.46,0.72,0.98,1.12];


figure(1)
hold on
plot(AR,CL)
scatter(AR,CL)
title('Simple Delta Wing Aspect Ratio vs C_L for \alpha = 20','interpreter','tex')
xlabel('AR')
ylabel('CL')
grid on