function [] = drawSim3(mainControlPointsY,mainControlPointsX,mainQuarterC,mainPanelGeomY,mainPanelGeomX,alpha,beta,wakeT,xw,yw,zw,xL,yL,zL,time,bPanels,cPanels,MAC,mode,anim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rootEnd = mainPanelGeomX(cPanels+1,bPanels+1);
rootStart = mainPanelGeomX(1,bPanels+1);

tipStart = mainPanelGeomX(1,1);
tipEnd = mainPanelGeomX(cPanels+1,1);


figure();
hold on
% view([0 alpha*(180/pi)]) 
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
    for j = 1:length(mainPanelGeomX(:,bPanels+1))
        plot([mainPanelGeomY(k) mainPanelGeomY(k)],[mainPanelGeomX(j,k) mainPanelGeomX(end,k)],'b')
        plot([mainPanelGeomY(k+1) mainPanelGeomY(k+1)],[mainPanelGeomX(j,k+1) mainPanelGeomX(end,k+1)],'b')
%         if j == 1
%             plot([mainPanelGeomY(k) mainPanelGeomY(k+1)],[mainPanelGeomX(j,k) mainPanelGeomX(j,k+1)],'b');            
%         else
            plot([mainPanelGeomY(k) mainPanelGeomY(k+1)],[mainPanelGeomX(j,k) mainPanelGeomX(j,k+1)],'b');  
%         end
    end
end

if mode == "fixedWake"
    for k = 1:length(mainPanelGeomY)
        plot3([mainPanelGeomY(k) mainPanelGeomY(k)+wakeT*tan(beta)],[mainPanelGeomX(end,k) -wakeT],[0 wakeT*tan(alpha)],'b')
    end
elseif mode == "freeWake"
       for k = 1:length(mainPanelGeomY)
           plot3(yw(:,k),xw(:,k),zw(:,k),'b')
           plot3(yL(:,k),xL(:,k),zL(:,k),'m');
%            plot3([yw(end,k) yw(end,k)+wakeT*tan(beta)],[xw(end,k) -wakeT],[zw(end,k) zw(end,k)+wakeT*tan(alpha)],'b')
       end    
end

panelCounter =1;
for g = 1:length(mainControlPointsY)
    for j = 1:cPanels
%         scatter(mainControlPointsY(g),mainControlPointsX(j,g),'g','x');
%             txt = num2str(panelCounter);
%              text(mainControlPointsY(g),mainQuarterC(j,g),txt)
        panelCounter = panelCounter +1;
    end
end

title('Wing and Wake Geometry');
xlabel('y');
ylabel('-x');
axis equal




if anim == "animate"  
   figure()
   hold on
   xlim([-40 40])
   ylim([-20 10])
   zlim([0 10])
   view([-37.5 30])
   filename = 'wakeAnimation.gif';
   for l = 1:length(time)+1
       plot3([0 0 mainPanelGeomY(1) mainPanelGeomY(1) 0],[rootStart rootEnd tipEnd tipStart rootStart],zeros(1,5),'r');
       plot3([0 mainPanelGeomY(end) mainPanelGeomY(end) 0 0],[rootStart tipStart tipEnd rootEnd rootStart],zeros(1,5),'r');  
       for k = 1:length(mainPanelGeomY)
            plot3(yw(1:l,k),xw(1:l,k),zw(1:l,k),'b')
       end
       drawnow
       M(l) = getframe;
       [imind,cm] = rgb2ind(frame2im(getframe),256);
       
       if l == 1
           imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0.25)
       else
           imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.25)
       end
       
   end  
   
   movie(M,1,4)
end
% ylim([-1*MAC 0.5])
% zlim([-0.1 0.2])
end

