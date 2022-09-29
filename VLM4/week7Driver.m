%Simple Rectangular
sections = 1;
taper = 1;
S = 100;
ARfix = 10;
b = sqrt(ARfix*S); %Span
rootC = S/b; %Root Chord
alphaIncRoot = 0*(pi/180); %Incidence at Root converted to radians
alphaIncTip = 0*(pi/180); %Incidence at Tip converted to radians
leLambda = linspace(0,60,7);
bPanels = 16;
cPanels = 16;
LE = 0;
rho =1;
camberLine = 0;
alpha = 5*(pi/180);
Uinf = 1;
for u = 1:length(leLambda)
    data = [taper;rootC;ARfix;leLambda(u)];
    [CLpot(u),~,~,CDiPot(u),~,~,CDiTreffz(u)] = vlm4aSweep(sections,data,bPanels,cPanels,alpha,Uinf,rho,LE);
end



figure()
plot(leLambda,CLpot);
title('C_L vs \Lambda','interpreter','tex')
xlabel('\Lambda');
ylabel('C_L');
grid on
saveas(gcf,'Fig1Rec.png')

figure()
hold on
plot(leLambda,CDiPot);
plot(leLambda,CDiTreffz);
title('C_Di vs \Lambda','interpreter','tex')
xlabel('\Lambda');
ylabel('C_Di');
legend('Wing Downwash','Treffz Plane Analysis');
grid on
saveas(gcf,'Fig2Rec.png')


%% Delta Wing

%Simple Delta
sections =1;
taper = 0; %Taper Ratio
S = 100; % Reference Area
ARfix = [0.25,0.5,1,1.5,2,3]; %Aspect Ratio
alphaIncRoot = 0*(pi/180); %Incidence at Root converted to radians
alphaIncTip = 0*(pi/180); %Incidence at Tip converted to radians
bPanels = 16;
cPanels =16;
alpha = 5*(pi/180);
Uinf = 1;
rho =1;
LE=0;
for u = 1:length(ARfix)
    b = sqrt(ARfix(u)*S); %Span
    rootC = S*2/b; %Root Chord
    leLambda= (90-atand((b/2)/rootC)); %Sweep angle 
    data = [taper;rootC;ARfix(u);leLambda];
    [CLpotDelta(u),~,~,CDiPotDelta(u),~,~,CDiTreffzDelta(u)] = vlm4aSweep(sections,data,bPanels,cPanels,alpha,Uinf,rho,LE);
end

figure()
plot(ARfix,CLpotDelta);
title('Simple Delta C_L vs AR','interpreter','tex')
xlabel('AR');
ylabel('C_L');
grid on
saveas(gcf,'Fig3Delta.png')

figure()
hold on
plot(ARfix,CDiPotDelta);
plot(ARfix,CDiTreffzDelta);
title('Simple Delta C_Di vs AR','interpreter','tex')
xlabel('AR');
ylabel('C_Di');
legend('Wing Downwash','Treffz Plane Analysis','Location','best');
grid on
saveas(gcf,'Fig4Delta.png')