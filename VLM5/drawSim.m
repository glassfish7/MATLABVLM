function [] = drawSim(mainControlPointsY,mainControlPointsX,mainQuarterC,mainPanelGeomY,mainPanelGeomX,alpha,beta,wakeT,xw,yw,zw,bPanels,cPanels,mode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rootEnd = mainPanelGeomX(cPanels+1,bPanels+1);
rootStart = mainPanelGeomX(1,bPanels+1);

tipStart = mainPanelGeomX(1,1);
tipEnd = mainPanelGeomX(cPanels+1,1);



figure()
hold on
       
plot([0 0 mainPanelGeomY(1) mainPanelGeomY(1) 0],[rootStart rootEnd tipEnd tipStart rootStart]);
plot([0 mainPanelGeomY(end) mainPanelGeomY(end) 0 0],[rootStart tipStart tipEnd rootEnd rootStart]);
        
for j = 1:length(mainPanelGeomX(:,bPanels+1))
   plot([0 mainPanelGeomY(1)],[mainPanelGeomX(j,bPanels+1) mainPanelGeomX(j,1)],'r');
   plot([0 mainPanelGeomY(end)],[mainPanelGeomX(j,bPanels+1) mainPanelGeomX(j,1)],'r');
end
        
for k = 1:length(mainPanelGeomY)
    plot([mainPanelGeomY(k) mainPanelGeomY(k)],[mainPanelGeomX(1,k) mainPanelGeomX(end,k)],'r');
end

for k = 1:length(mainControlPointsY)
    for j = 1:length(mainQuarterC(:,bPanels+1))
        plot([mainPanelGeomY(k) mainPanelGeomY(k)],[mainQuarterC(j,k) mainPanelGeomX(end,k)],'m')
        plot([mainPanelGeomY(k+1) mainPanelGeomY(k+1)],[mainQuarterC(j,k+1) mainPanelGeomX(end,k+1)],'m')
        plot([mainPanelGeomY(k) mainPanelGeomY(k+1)],[mainQuarterC(j,k) mainQuarterC(j,k+1)],'b');     
    end
end

if mode == "freeWake"
    for k = 1:length(mainPanelGeomY)
        plot3([mainPanelGeomY(k) mainPanelGeomY(k)+wakeT*tan(beta)],[mainPanelGeomX(end,k) -wakeT],[0 wakeT*tan(alpha)],'b')
    end
elseif mode == "realWake"
       for k = 1:length(mainPanelGeomY)
           plot3(yw(:,k),xw(:,k),zw(:,k),'b')
       end
    
end


panelCounter =1;
for g = 1:length(mainControlPointsY)
    for j = 1:cPanels
        scatter(mainControlPointsY(g),mainControlPointsX(j,g),'g','x');
%             txt = num2str(panelCounter);
%              text(mainControlPointsY(g),mainQuarterC(j,g),txt)
        panelCounter = panelCounter +1;
    end
end

title('Wing and Wake Geometry');
xlabel('y');
ylabel('-x');
axis equal


end

