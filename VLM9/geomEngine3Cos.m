function [mainControlPointsY,mainControlPointsX,mainQuarterC,maintquarterPointsX,mainPanelGeomY,mainPanelGeomX,mainVControlPointsX,Stot,span,Npanels,K,CAve,chordDistGeom,chordDistCont] = geomEngine3(sections,data,bPanelsL,cPanels,plotCont)
%geomEngine2 Generates Panels and Control Points for VLM
%   Detailed explanation goes here

if size(data) ~= [4,sections]
    error('Data mismatch')
end

if length(bPanelsL) ~= sections
    error('Define Panels')
end


tquartX = cell(1,sections);
panelQC = cell(1,sections);
panelGY = cell(1,sections);
panelGX = cell(1,sections);


Stot = 0;
span =0;
Npanels = 2*sum(bPanelsL)*cPanels;

for sec = 1:sections
 
%Select panel number set 
bPanels = bPanelsL(sec);


%GeomEngine Baseline

if sec ==1
    taper = data(1,sec); %Taper Ratio
    rootC = data(2,sec); % Root Chord
    AR = data(3,sec); %Aspect Ratio
    tipC = rootC*taper; % Tip Chord
    S =  AR*((tipC+rootC)/2)^2; %Area
    b = sqrt(AR*S); %Span
    leLambda = data(4,sec)*(pi/180);%LE Sweep angle in degrees/converted to Radians
    Stot = Stot +S;
    span = span + b;
    mainrootC = rootC;
    maintipC = tipC;
else
    taper = data(1,sec); %Taper Ratio
    rootC = abs(panelGX{sec-1}(1,1) - panelGX{sec-1}(cPanels+1,1));
    AR = data(3,sec); %Aspect Ratio
    tipC = rootC*taper; % Tip Chord
    S =  AR*((tipC+rootC)/2)^2; %Area
    b = sqrt(AR*S); %Span
    leLambda = data(4,sec)*(pi/180);%LE Sweep angle in degrees/converted to Radians
    Stot = Stot +S;
    span = span + b;
    if sec == sections
        maintipC = tipC;
    end
    
end    


CAve = (maintipC + mainrootC)/2;


%Calculate Root and Tip ending and starting points
if sec ==1
    rootEnd = 0;
    rootStart = rootC + rootEnd;

    tipStart = rootStart - (b/2)*tan(leLambda);
    tipEnd = tipStart - tipC;
else
   rootEnd =  panelGX{sec-1}(cPanels+1,1);
   rootStart =   rootC + rootEnd;
   
   tipStart = rootStart - (b/2)*tan(leLambda);
   tipEnd = tipStart - tipC;
end



%Generate Panels

%Descretize wing root and wing tip into N+1 points  
% rootChordPoints = rootStart:-(rootStart-rootEnd)/cPanels:rootEnd;
% tipChordPoints =  tipStart:-(tipStart-tipEnd)/cPanels:tipEnd;

rootChordPoints = zeros(1,cPanels+1);
tipChordPoints = zeros(1,cPanels+1);

for space = 0:cPanels
    rootChordPoints(space+1) = rootEnd + (rootStart-rootEnd)*cos((pi/2)*(space/cPanels));
    tipChordPoints(space+1) = tipEnd + (tipStart-tipEnd)*cos((pi/2)*(space/cPanels));
end

if isempty(tipChordPoints)
    tipChordPoints = tipStart*ones(1,length(rootChordPoints));
end


%Divide wing spanwise into 2N+1 stations
K(sec) = b/(2*bPanels);
if sec == 1
    panelGeomY = -b/2:K(sec):b/2;
else
    panelGeomY = zeros(1,2*bPanels);
    panelGeomY(1:bPanels) = (panelGY{1,sec-1}(1)-(b/2)):K(sec):(panelGY{1,sec-1}(1)-K(sec));
    panelGeomY(bPanels+1:2*bPanels) = (panelGY{1,sec-1}(end)+K(sec)):K(sec):(panelGY{1,sec-1}(end)+(b/2));
end


    %Find wing "centerline" index
    if sec ==1
        centreIndex = bPanels +1;
    else 
        centreIndex = bPanels;
    end
    
    %Create panel x coordiantes
    panelGeomX  = zeros(length(rootChordPoints),length(panelGeomY));
    if sec == 1  
        for i = 1:length(rootChordPoints)
            %Left Wing
            panelGeomX(i,1:centreIndex) = rootChordPoints(i) + ((tipChordPoints(i)-rootChordPoints(i))/(-b/2))*(panelGeomY(1:centreIndex));
            %Right Wing
            panelGeomX(i,centreIndex+1:end) = rootChordPoints(i) + ((tipChordPoints(i)-rootChordPoints(i))/(b/2))*(panelGeomY(centreIndex+1:end));

        end
    else
        for i = 1:length(rootChordPoints)
            %Left Wing
            panelGeomX(i,1:centreIndex) = rootChordPoints(i) + ((tipChordPoints(i)-rootChordPoints(i))/(-b/2))*((panelGeomY(1:centreIndex))-panelGY{1,sec-1}(1));
            %Right Wing
            panelGeomX(i,centreIndex+1:end) = rootChordPoints(i) + ((tipChordPoints(i)-rootChordPoints(i))/(b/2))*(panelGeomY(centreIndex+1:end) - panelGY{1,sec-1}(end));

        end
    end
    

    %Create quarter chord points, three quarter chord points for each panel
    panelQuarterC = zeros(cPanels,length(panelGeomY));
    tquarterPointsX = zeros(cPanels,length(panelGeomY));
    for i = 1:length(panelGeomY)
        for j = 1:cPanels
            panelQuarterC(j,i) = panelGeomX(j,i) + (panelGeomX(j+1,i)-panelGeomX(j,i))/4;
            tquarterPointsX(j,i) = panelGeomX(j,i) + 3*(panelGeomX(j+1,i)-panelGeomX(j,i))/4;
        end
    end
       

panelQC(1,sec) = {panelQuarterC};
panelGY(1,sec) = {panelGeomY};
tquartX(1,sec) = {tquarterPointsX};
panelGX(1,sec) = {panelGeomX};


end

%Stitch the results
for i = 1:sections
    if i == 1
        mainQuarterC = panelQC{1,i};
        mainPanelGeomY = panelGY{1,i};
        mainPanelGeomX = panelGX{1,i};
        maintquarterPointsX = tquartX{1,i};
    else
        mainQuarterC = [panelQC{1,i}(:,1:bPanelsL(i)) mainQuarterC panelQC{1,i}(:,bPanelsL(i)+1:end)];
        maintquarterPointsX = [tquartX{1,i}(:,1:bPanelsL(i)) maintquarterPointsX tquartX{1,i}(:,bPanelsL(i)+1:end)];
        mainPanelGeomY = [panelGY{1,i}(:,1:bPanelsL(i)) mainPanelGeomY panelGY{1,i}(:,bPanelsL(i)+1:end)];
        mainPanelGeomX = [panelGX{1,i}(:,1:bPanelsL(i)) mainPanelGeomX panelGX{1,i}(:,bPanelsL(i)+1:end)];
    end
    
end

%Use sitched geometry to calculate control point locations
%and chord distribution
mainControlPointsX = zeros(cPanels,2*sum(bPanelsL));
mainVControlPointsX = zeros(cPanels,2*sum(bPanelsL));
mainControlPointsY= zeros(1,2*sum(bPanelsL));
chordDistGeom = zeros(1,length(mainPanelGeomY));
chordDistCont = zeros(1,length(mainControlPointsY));

for i = 1:length(mainPanelGeomY)
    chordDistGeom(i)  = mainPanelGeomX(1,i)-mainPanelGeomX(end,i);
end

for i= 1:length(mainPanelGeomY)-1
   for j = 1:cPanels
         mainControlPointsX(j,i) = (maintquarterPointsX(j,i+1)+maintquarterPointsX(j,i))/2;
         mainVControlPointsX(j,i) = (mainQuarterC(j,i+1)+mainQuarterC(j,i))/2;
   end
     mainControlPointsY(i) = mainPanelGeomY(i) + (mainPanelGeomY(i+1)-mainPanelGeomY(i))/2;
     chordDistCont(i) = (chordDistGeom(i+1) + chordDistGeom(i))/2;
end
      

%Plotting Section

if plotCont ~= "NoPlot"
%Plot Planform
figure()
hold on
if sec ==1
    plot([0 0 -b/2 -b/2 0 0 b/2 b/2 0],[rootStart rootEnd tipEnd tipStart rootStart rootStart tipStart tipEnd rootEnd]); 
else
    plot([panelGY{1,sec-1}(1) panelGY{1,sec-1}(1) panelGY{1,sec-1}(1)-b/2 panelGY{1,sec-1}(1)-b/2 panelGY{1,sec-1}(1) panelGY{1,sec-1}(end) panelGY{1,sec-1}(end)+b/2 panelGY{1,sec-1}(end)+b/2 panelGY{1,sec-1}(end)],[rootStart rootEnd tipEnd tipStart rootStart rootStart tipStart tipEnd rootEnd]);
end
title('Section Dimensions and Shape')
xlabel('y');
ylabel('-x');
axis equal


figure()
hold on
for i = 1:sections
    
    if i == 1
        rootEnd = panelGX{i}(cPanels+1,bPanelsL(i)+1);
        rootStart = panelGX{i}(1,bPanelsL(i)+1);

        tipStart = panelGX{i}(1,1);
        tipEnd = panelGX{i}(cPanels+1,1);
        
        plot([0 0 panelGY{1,i}(1) panelGY{1,i}(1) 0],[rootStart rootEnd tipEnd tipStart rootStart]);
        plot([0 panelGY{1,i}(end) panelGY{1,i}(end) 0 0],[rootStart tipStart tipEnd rootEnd rootStart]);
        
        for j = 1:length(panelGX{i}(:,bPanelsL(i)+1))
            plot([0 panelGY{1,i}(1)],[panelGX{i}(j,bPanelsL(i)+1) panelGX{i}(j,1)],'r');
            plot([0 panelGY{1,i}(end)],[panelGX{i}(j,bPanelsL(i)+1) panelGX{i}(j,1)],'r');
        end
        
        for k = 1:length(panelGY{1,i})
            plot([panelGY{1,i}(k) panelGY{1,i}(k)],[panelGX{i}(1,k) panelGX{i}(end,k)],'r');
        end

    else
        rootEnd = panelGX{i-1}(cPanels+1,1);
        rootStart = panelGX{i-1}(1,1);

        tipStart = panelGX{i}(1,1);
        tipEnd = panelGX{i}(cPanels+1,1);
        
        plot([panelGY{1,i-1}(1) panelGY{1,i-1}(1) panelGY{1,i}(1) panelGY{1,i}(1) panelGY{1,i-1}(1)],[rootStart rootEnd tipEnd tipStart rootStart]);
        plot([panelGY{1,i-1}(end) panelGY{1,i}(end) panelGY{1,i}(end) panelGY{1,i-1}(end) panelGY{1,i-1}(end)],[rootStart tipStart tipEnd rootEnd rootStart]);   
        
        for j = 1:length(panelGX{i-1}(:,1))
            plot([panelGY{1,i-1}(1) panelGY{1,i}(1)],[panelGX{i-1}(j,1) panelGX{i}(j,1)],'r');
            plot([panelGY{1,i-1}(end) panelGY{1,i}(end)],[panelGX{i-1}(j,1) panelGX{i}(j,1)],'r');
        end
        
        for k = 1:length(panelGY{1,i})
            plot([panelGY{1,i}(k) panelGY{1,i}(k)],[panelGX{i}(1,k) panelGX{i}(end,k)],'r');
        end
    end
    
    panelCounter =1;
    for g = 1:length(mainControlPointsY)
        for j = 1:cPanels
            scatter(mainControlPointsY(g),mainControlPointsX(j,g),'b','x');
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

end

