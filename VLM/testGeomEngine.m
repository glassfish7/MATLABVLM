sections = 3;
data = [0.75,0.5,0.1;1,1,1;5,4,3;0,0,0];
bPanels =4;
cPanels =4;




if size(data) ~= [4,sections]
    error('Data mismatch')
end


contPY = cell(1,sections);
contPX = cell(1,sections);
panelQC = cell(1,sections);
panelGY = cell(1,sections);
contPVX = cell(1,sections);
panelGX = cell(1,sections);


for sec = 1:sections
 

%GeomEngine Baseline

if sec ==1
    taper = data(1,sec); %Taper Ratio
    rootC = data(2,sec); % Root Chord
    AR = data(3,sec); %Aspect Ratio
    tipC = rootC*taper; % Tip Chord
    S =  AR*((tipC+rootC)/2)^2; %Area
    b = sqrt(AR*S); %Span
    %leLambda= (90-atand((b/2)/rootC))*(pi/180); 
    leLambda = data(4,sec)*(pi/180);%LE Sweep angle in degrees/converted to Radians
else
    taper = data(1,sec); %Taper Ratio
    rootC = abs(panelGX{sec-1}(1,1) - panelGX{sec-1}(cPanels+1,1));

    AR = data(3,sec); %Aspect Ratio
    tipC = rootC*taper; % Tip Chord
    S =  AR*((tipC+rootC)/2)^2; %Area
    b = sqrt(AR*S); %Span
    %leLambda= (90-atand((b/2)/rootC))*(pi/180); 
    leLambda = data(4,sec)*(pi/180);%LE Sweep angle in degrees/converted to Radians

end    
%Calculate Root and Tip ending and starting points
if sec ==1
    rootEnd = 0;
    rootStart = rootC + rootEnd;

    tipStart = rootStart - (b/2)*tan(leLambda);
    tipEnd = tipStart - tipC;
else
   rootEnd =  panelGX{sec-1}(cPanels+1,1);
   disp(rootEnd)
   rootStart =   rootC + rootEnd;
   
   tipStart = rootStart - (b/2)*tan(leLambda);
   tipEnd = tipStart - tipC;
end

%Plot Planform
figure(1)
hold on
if sec ==1
    plot([0 0 -b/2 -b/2 0 0 b/2 b/2 0],[rootStart rootEnd tipEnd tipStart rootStart rootStart tipStart tipEnd rootEnd]); 
else
    plot([panelGY{1,sec-1}(1) panelGY{1,sec-1}(1) panelGY{1,sec-1}(1)-b/2 panelGY{1,sec-1}(1)-b/2 panelGY{1,sec-1}(1) panelGY{1,sec-1}(end) panelGY{1,sec-1}(end)+b/2 panelGY{1,sec-1}(end)+b/2 panelGY{1,sec-1}(end)],[rootStart rootEnd tipEnd tipStart rootStart rootStart tipStart tipEnd rootEnd]);
end
title('Section Dimentions and Shape')
xlabel('y');
ylabel('-x');
axis equal



%Generate Panels

%Descretize wing root and wing tip into N+1 points  
rootChordPoints = rootStart:-(rootStart-rootEnd)/cPanels:rootEnd;
tipChordPoints =  tipStart:-(tipStart-tipEnd)/cPanels:tipEnd;

if isempty(tipChordPoints)
    tipChordPoints = tipStart*ones(1,length(rootChordPoints));
end


%Divide wing spanwise into 2N+1 stations
K = b/(2*bPanels);
if sec == 1
    panelGeomY = -b/2:K:b/2;
else
    panelGeomY = zeros(1,2*bPanels);
    panelGeomY(1:bPanels) = (panelGY{1,sec-1}(1)-(b/2)):K:(panelGY{1,sec-1}(1)-K);
    panelGeomY(bPanels+1:2*bPanels) = (panelGY{1,sec-1}(end)+K):K:(panelGY{1,sec-1}(end)+(b/2));
end


    %Find wing "centerline" index
    if sec ==1
        centreIndex = bPanels +1;
    else 
        centreIndex = bPanels;
    end
    
    if sec == 1
        %Create panel x coordiantes
        panelGeomX  = zeros(length(rootChordPoints),length(panelGeomY));
        for i = 1:length(rootChordPoints)
            %Left Wing
            panelGeomX(i,1:centreIndex) = rootChordPoints(i) + ((tipChordPoints(i)-rootChordPoints(i))/(-b/2))*(panelGeomY(1:centreIndex));
            %Right Wing
            panelGeomX(i,centreIndex+1:end) = rootChordPoints(i) + ((tipChordPoints(i)-rootChordPoints(i))/(b/2))*(panelGeomY(centreIndex+1:end));

        end
    else
        panelGeomX  = zeros(length(rootChordPoints),length(panelGeomY));
        for i = 1:length(rootChordPoints)
            %Left Wing
            panelGeomX(i,1:centreIndex) = rootChordPoints(i) + ((tipChordPoints(i)-rootChordPoints(i))/(-b/2))*((panelGeomY(1:centreIndex))-panelGY{1,sec-1}(1));
            %Right Wing
            panelGeomX(i,centreIndex+1:end) = rootChordPoints(i) + ((tipChordPoints(i)-rootChordPoints(i))/(b/2))*(panelGeomY(centreIndex+1:end) - panelGY{1,sec-1}(end));

        end
    end
    

    %Create quarter chord points for each panel
    panelQuarterC = zeros(cPanels,length(panelGeomY));
    for i = 1:length(panelGeomY)
        for j = 1:cPanels
            panelQuarterC(j,i) = panelGeomX(j,i) + (panelGeomX(j+1,i)-panelGeomX(j,i))/4;
        end
    end


    %Create three quarter chord points for each panel
    tquarterPointsX = zeros(cPanels,length(panelGeomY));
    for i = 1:length(panelGeomY)
        for j = 1:cPanels
            tquarterPointsX(j,i) = panelGeomX(j,i) + 3*(panelGeomX(j+1,i)-panelGeomX(j,i))/4;
        end
    end


    %Create control point x coordinate for each panel
    
    if sec ==1
        controlPointsX = zeros(cPanels,length(panelGeomY)-1);
        for i= 1:length(panelGeomY)-1
            for j = 1:cPanels
                controlPointsX(j,i) = (tquarterPointsX(j,i+1)+tquarterPointsX(j,i))/2;
            end
        end
    else
        controlPointsX = zeros(cPanels,length(panelGeomY)-1);
        for i= 1:length(panelGeomY)-1
            for j = 1:cPanels
                if i == bPanels
                    controlPointsX(j,i) = (tquartX{1,sec-1}(j,1)+tquarterPointsX(j,i))/2;
                elseif i == bPanels +1
                    controlPointsX(j,i) = (tquarterPointsX(j,i+1)+tquartX{1,sec-1}(j,end))/2;
                else
                    controlPointsX(j,i) = (tquarterPointsX(j,i+1)+tquarterPointsX(j,i))/2;
                end
            end
        end
        
    end
    
    if sec == 1
        %Create on vortex control points for vortex lift calculations
        VcontrolPointsX = zeros(cPanels,length(panelGeomY)-1);
        for i= 1:length(panelGeomY)-1
            for j = 1:cPanels
                VcontrolPointsX(j,i) = (panelQuarterC(j,i+1)+panelQuarterC(j,i))/2;
            end
        end
    else
        VcontrolPointsX = zeros(cPanels,length(panelGeomY)-1);
        for i= 1:length(panelGeomY)-1
            for j = 1:cPanels
                
                if i == bPanels
                    VcontrolPointsX(j,i) = (panelQC{1,sec-1}(j,1)+panelQuarterC(j,i))/2;
                elseif i == bPanels+1
                    VcontrolPointsX(j,i) = (panelQuarterC(j,i+1)+panelQC{1,sec-1}(j,end))/2;
                else
                    VcontrolPointsX(j,i) = (panelQuarterC(j,i+1)+panelQuarterC(j,i))/2;
                end
            end
        end
        
    end

    %Create control point y coordinate for each panel
    
    if sec == 1
        controlPointsY= zeros(1,length(panelGeomY)-1);
        for i = 1:length(panelGeomY)-1
            controlPointsY(i) = panelGeomY(i) + (panelGeomY(i+1)-panelGeomY(i))/2;
        end
    else
        
        controlPointsY= zeros(1,length(panelGeomY)-1);
        for i = 1:length(panelGeomY)-1
            if i <= bPanels
                controlPointsY(i) = panelGeomY(i) + K/2;
            else
                controlPointsY(i) = panelGeomY(i) - K/2;
            end
        end
        
        
        
    end
    

    
contPY(1,sec) = {controlPointsY};
contPX(1,sec) = {controlPointsX};
panelQC(1,sec) = {panelQuarterC};
panelGY(1,sec) = {panelGeomY};
tquartX(1,sec) = {tquarterPointsX};
panelGX(1,sec) = {panelGeomX};
contPVX(1,sec) = {VcontrolPointsX};


end

%Stitch the results
% mainControlPointsY = zeros(1,sections*2*bPanels);
% mainControlPointsX = zeros(cpanels,sections*2*bPanels);
% mainQuarterC = zeros(cpanels,sections*2*bPanels);
% mainPanelGeomY = zeros(1,sections*2*bPanels+1);
% mainPanelGeomX = zeros(cpanels+1,sections*2*bPanels+1);
% mainVControlPointsX = zeros(cpanels,sections*2*bPanels);
for i = 1:sections
    if i == 1
        mainControlPointsY = contPY{1,i};
        mainControlPointsX = contPX{1,i};
        mainQuarterC = panelQC{1,i};
        mainPanelGeomY = panelGY{1,i};
        mainPanelGeomX = panelGX{1,i};
        mainVControlPointsX = contPVX{1,i};
    else
        mainControlPointsY = [contPY{1,i}(:,1:bPanels) mainControlPointsY contPY{1,i}(:,bPanels+1:end)];
        mainControlPointsX = [contPX{1,i}(:,1:bPanels) mainControlPointsX contPX{1,i}(:,bPanels+1:end)];
        mainQuarterC = [panelQC{1,i}(:,1:bPanels) mainQuarterC panelQC{1,i}(:,bPanels+1:end)];
        mainPanelGeomY = [panelGY{1,i}(:,1:bPanels) mainPanelGeomY panelGY{1,i}(:,bPanels+1:end)];
        mainPanelGeomX = [panelGX{1,i}(:,1:bPanels) mainPanelGeomX panelGX{1,i}(:,bPanels+1:end)];
        mainVControlPointsX = [contPVX{1,i}(:,1:bPanels) mainVControlPointsX contPVX{1,i}(:,bPanels+1:end)];
    end
    
end