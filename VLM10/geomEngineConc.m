function [controlPointsY,controlPointsX,panelQuarterC,panelQuarterCX,panelHalfCY,tquarterPointsX,panelGeomY,panelGeomX,vControlPointsX,VcontrolPointsY,Stot,span,Npanels,K,CAve,chordDistGeom,chordDistCont] = geomEngineConcSpec(~,data,bPanelsL,cPanels,plotCont)
%geomEngine3Conc Generates Concial Panels and Control Points for VLM
%   Detailed explanation goes here

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



%Generate Panels

%Descretize wing leading edge and wing trailing edge into N+1 points  
% trailingEdgePoints = linspace(-b/2,b/2,2*bPanels+1);
% leadingEdgePoints = zeros(1,length(trailingEdgePoints));


%Divide wing chordwise into 2N+1 stations
K = rootC/(cPanels);
panelGeomX = (rootStart:-K:rootEnd)';


    
%Create panel y coordiantes
panelGeomY  = zeros(length(panelGeomX),length(trailingEdgePoints));
for i = 1:length(trailingEdgePoints)
        %Left Wing
        panelGeomY(:,i) = trailingEdgePoints(i) + ((leadingEdgePoints(i)-trailingEdgePoints(i))/(rootC))*panelGeomX;
end

    

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
tGeomY  = zeros(cPanels,length(trailingEdgePoints));
qGeomY = zeros(cPanels,length(trailingEdgePoints));
panelHalfCY = zeros(cPanels,length(trailingEdgePoints));
panelQuarterC = zeros(cPanels,2*bPanels+1);
for i = 1:length(trailingEdgePoints)
        %Left Wing
        tGeomY(:,i) = trailingEdgePoints(i) + ((leadingEdgePoints(i)-trailingEdgePoints(i))/(rootC))*tquarterPointsX(:,1);
        qGeomY(:,i) = trailingEdgePoints(i) + ((leadingEdgePoints(i)-trailingEdgePoints(i))/(rootC))*qquarterPointsX(:,1);
        panelQuarterC(:,i) = trailingEdgePoints(i) + ((leadingEdgePoints(i)-trailingEdgePoints(i))/(rootC))*panelQuarterCX(:,1);
        panelHalfCY(:,i) = trailingEdgePoints(i) + ((leadingEdgePoints(i)-trailingEdgePoints(i))/(rootC))*halfPointsX(:,1);
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

    
    plot([0 0 panelGeomY(end,1) panelGeomY(end,1) 0],[rootStart rootEnd tipEnd tipStart rootStart]);
    plot([0 panelGeomY(end,end) panelGeomY(end,end) 0 0],[rootStart tipStart tipEnd rootEnd rootStart]);

    for j = 1:length(trailingEdgePoints)
        plot([trailingEdgePoints(j) 0],[panelGeomX(end) panelGeomX(1)],'r')
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

