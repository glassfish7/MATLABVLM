function [controlPointsY,controlPointsX,panelQuarterC,panelQuarterCX,panelHalfCY,tquarterPointsX,panelGeomY,panelGeomX,vControlPointsX,VcontrolPointsY,Stot,span,Npanels,K,CAve,chordDistGeom,chordDistCont] = geomEngineConcSpec(~,data,bPanelsL,cPanels,plotCont)
%geomEngine3Conc Generates Concial Panels and Control Points for VLM
%   Detailed explanation goes here
K =1;
chordDistGeom = 0;
chordDistCont = 0;
%GeomEngine Baseline

    taper = data(1); %Taper Ratio
    rootC = data(2); % Root Chord
    AR = data(3); %Aspect Ratio
    tipC = rootC*taper; % Tip Chord
    S =  AR*((tipC+rootC)/2)^2; %Area
    b = sqrt(AR*S); %Span
    leLambda = data(4)*(pi/180);%LE Sweep angle in degrees/converted to Radians
    Stot = S;
    span = b;
    mainrootC = rootC;
    maintipC = tipC;
    bPanels = bPanelsL;
    Npanels = cPanels*2*bPanels;
    
CAve = (maintipC + mainrootC)/2;


%Calculate Root and Tip ending and starting points
rootEnd = 0;
rootStart = rootC + rootEnd;

tipStart = rootStart - (b/2)*tan(leLambda);
tipEnd = tipStart - tipC;

cPanels = 24;
bPanels = 12;
%Generate Panels
%Hardcoded
K2 = 5.68294/(12);
K1 = 8.76706/(12);
panelGeomX = [14.45:-K1:5.68294 5.68294-K2:-K2:0]';

slope1 = 2.25/5.68294;
slope2 = 2.75/(14.45-5.68294);

levels1 = -2.75 + slope1*([5.68294-K2:-K2:0]-5.68294)';
levels2 = slope2*([14.45:-K1:5.68294]-14.45)';
levels  = [levels2;levels1];

points = zeros(25,25);

for i = 1:length(levels)
    panelGeomY(i,:) = linspace(-levels(i),levels(i),25);
end

%Descretize wing leading edge and wing trailing edge into N+1 points  
% trailingEdgePoints = linspace(-b/2,b/2,2*bPanels+1);
% leadingEdgePoints = zeros(1,length(trailingEdgePoints));



    

%Create quarter chord points, three quarter chord points for each panel
panelQuarterCX = zeros(cPanels,2*bPanels);
tquarterPointsX = zeros(cPanels,2*bPanels);
qquarterPointsX = zeros(cPanels,2*bPanels);
halfPointsX = zeros(cPanels,2*bPanels);
for i = 1:2*bPanels
    for j = 1:cPanels
        panelQuarterCX(j,i) = panelGeomX(j) + (panelGeomX(j+1)-panelGeomX(j))/4;
        tquarterPointsX(j,i) = panelGeomX(j) + 3*(panelGeomX(j+1)-panelGeomX(j))/4;
        qquarterPointsX(j,i) = panelGeomX(j) + (panelGeomX(j+1)-panelGeomX(j))/4;
        halfPointsX(j,i) = panelGeomX(j) + (panelGeomX(j+1)-panelGeomX(j))/2;
    end
end
       
%Special Control Point Y Algorithm. Concical Panelling only
tGeomY  = zeros(cPanels,2*bPanels+1);
qGeomY = zeros(cPanels,2*bPanels+1);
panelHalfCY = zeros(cPanels,2*bPanels+1);
panelQuarterC = zeros(cPanels,2*bPanels+1);
for i = 1:2*bPanels+1
    for j = 1:12
        tGeomY(j,i) = panelGeomY(j,i) + ((panelGeomY(j+1,i)-panelGeomY(j,i))/(K1))*(3*K1/4);
        qGeomY(j,i) = panelGeomY(j,i) + ((panelGeomY(j+1,i)-panelGeomY(j,i))/(K1))*(K1/4);
        panelQuarterC(j,i) = panelGeomY(j,i) + ((panelGeomY(j+1,i)-panelGeomY(j,i))/(K1))*(K1/4);
        panelHalfCY(j,i) = panelGeomY(j,i) + ((panelGeomY(j+1,i)-panelGeomY(j,i))/(K1))*(K1/2);
    end
end

for i = 1:2*bPanels+1
    for j = 13:24
        tGeomY(j,i) = panelGeomY(j,i) + ((panelGeomY(j+1,i)-panelGeomY(j,i))/(K2))*(3*K2/4);
        qGeomY(j,i) = panelGeomY(j,i) + ((panelGeomY(j+1,i)-panelGeomY(j,i))/(K2))*(K2/4);
        panelQuarterC(j,i) = panelGeomY(j,i) + ((panelGeomY(j+1,i)-panelGeomY(j,i))/(K2))*(K2/4);
        panelHalfCY(j,i) = panelGeomY(j,i) + ((panelGeomY(j+1,i)-panelGeomY(j,i))/(K2))*(K2/2);
    end
end




%Use sitched geometry to calculate control point locations
%and chord distribution
controlPointsX = tquarterPointsX;
vControlPointsX = qquarterPointsX;
controlPointsY= zeros(cPanels,2*sum(bPanels));
VcontrolPointsY = zeros(cPanels,2*sum(bPanels));
% for i = 1:length(mainPanelGeomY)
%     chordDistGeom(i)  = mainPanelGeomX(1,i)-mainPanelGeomX(end,i);
% end
for j = 1:cPanels
    for i= 1:2*bPanels
         controlPointsY(j,i) = tGeomY(j,i) + (tGeomY(j,i+1)-tGeomY(j,i))/2;
         VcontrolPointsY(j,i) = qGeomY(j,i) + (qGeomY(j,i+1)-qGeomY(j,i))/2;
    end
end
      

%Plotting Section

if plotCont ~= "NoPlot"

figure()
hold on

    
%     plot([0 0 panelGeomY(end,1) panelGeomY(end,1) 0],[rootStart rootEnd tipEnd tipStart rootStart]);
%     plot([0 panelGeomY(end,end) panelGeomY(end,end) 0 0],[rootStart tipStart tipEnd rootEnd rootStart]);

    for j = 1:cPanels
        for i = 1:2*bPanels+1
        plot([panelGeomY(j,i) panelGeomY(j+1,i)],[panelGeomX(j) panelGeomX(j+1)],'r');
        end
    end

    for k = 1:length(panelGeomX)
        plot([panelGeomY(k,1) panelGeomY(k,end)],[panelGeomX(k) panelGeomX(k)],'r');
    end

    
    panelCounter =1;
    for g = 1:2*bPanels
        for j = 1:cPanels
            scatter(controlPointsY(j,g),controlPointsX(j,g),'b','x');
%             txt = num2str(panelCounter);
%              text(mainControlPointsY(g),mainQuarterC(j,g),txt)
            panelCounter = panelCounter +1;
        end
    end
end
title('Descretized Wing Geometry Output');
xlabel('y');
ylabel('-x');
axis equal

end

