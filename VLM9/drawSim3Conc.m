function [] = drawSim3Conc(mainControlPointsY,mainControlPointsX,mainQuarterC,mainQuarterCX,mainPanelGeomY,mainPanelGeomX,alpha,beta,wakeT,xw,yw,zw,xL,yL,zL,time,bPanels,cPanels,MAC,mode,anim)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rootEnd = mainPanelGeomX(end);
rootStart = mainPanelGeomX(1);

tipStart = mainPanelGeomX(end);
tipEnd = mainPanelGeomX(end);


figure();
hold on
% view([0 0]) 
% plot([0 0 mainPanelGeomY(end,1) mainPanelGeomY(end,1) 0],[rootStart rootEnd tipEnd tipStart rootStart]);
% plot([0 mainPanelGeomY(end,end) mainPanelGeomY(end,end) 0 0],[rootStart tipStart tipEnd rootEnd rootStart]);
         
for j = 1:2*bPanels+1
    plot([mainPanelGeomY(end,j) 0],[mainPanelGeomX(end) mainPanelGeomX(1)],'r')
end

for k = 1:length(mainPanelGeomX)
    plot([mainPanelGeomY(k,1) mainPanelGeomY(k,end)],[mainPanelGeomX(k) mainPanelGeomX(k)],'r');
end

for k = 1:2*bPanels
    for j = 1:cPanels

        plot([mainPanelGeomY(end,k) mainQuarterC(j,k)],[mainPanelGeomX(end) mainQuarterCX(j,k)],'b')
        plot([mainPanelGeomY(end,k+1) mainQuarterC(j,k+1)],[mainPanelGeomX(end) mainQuarterCX(j,k)],'b')

    end
end

if mode == "fixedWake"
    for k = 1:length(mainPanelGeomY)
        plot3([mainPanelGeomY(k) mainPanelGeomY(k)+wakeT*tan(beta)],[mainPanelGeomX(end,k) -wakeT],[0 wakeT*tan(alpha)],'b')
    end
elseif mode == "freeWake"
       for k = 1:2*bPanels+1
           plot3(yw(:,k),xw(:,k),zw(:,k),'b')

% %            plot3([yw(end,k) yw(end,k)+wakeT*tan(beta)],[xw(end,k) -wakeT],[zw(end,k) zw(end,k)+wakeT*tan(alpha)],'b')
       end    
%             
          for k = 1:cPanels-1
              plot3(yL(:,k),xL(:,k),zL(:,k),'m');
          end
          for k = cPanels+2:2*cPanels

              plot3(yL(:,k),xL(:,k),zL(:,k),'m');
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

title('Wake Shape');
xlabel('y/c_r');
zlabel('-z/c_r');
% set ( gca, 'zdir', 'reverse' )
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
% ylim([-1*20 25])
% zlim([-1 3])
end

