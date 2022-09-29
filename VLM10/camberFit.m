function [camberDrep,camberMat,camber] = camberFit(camberLine,cPanels,bPanels,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%Set Airfoil sections
%Set normalized control point locations
normControlPoints = linspace(0,1,cPanels+1);
normControlPoints  = normControlPoints + (3/(4*cPanels));
normControlPoints  = normControlPoints(1:end-1);
%Assign camber line to each chordwise wing section
%Take derivative with central difference
%3 approaches: 4-digit NACA, custom camber line, and no camber
if isa(camberLine,'char')
    if length(camberLine) ~= 4
        error('Please enter a 4 digit NACA airfoil');
    else
        camber = zeros(1,length(normControlPoints));
        camberD = zeros(1,length(normControlPoints));
        m = str2double(camberLine(1))/100;
        p = str2double(camberLine(2))/10;
        for i = 1:length(normControlPoints)
            if normControlPoints <= p
                camber(i) = (m/p^2)*(2*p*normControlPoints(i) - normControlPoints(i).^2);
                camberD(i) = (2*m/p^2)*(p - normControlPoints(i));
            else
                camber(i) = (m/(1-p)^2)*(1-2*p+2*p*normControlPoints(i) - normControlPoints(i).^2);
                camberD(i) = (2*m/(1-p)^2)*(p - normControlPoints(i));
            end
        end
    end
    camberMat  = repmat(camberD',1,2*bPanels);
    camberDrep = repmat(camberD',2*bPanels,1);
elseif camberLine == 0
    camber = 0;
    camberDrep = zeros(N,1);
    camberMat = zeros(cPanels,2*bPanels);
else 
    %Custom camber line functionality. 
    %Multi section airfoil functionality to be added later
    [camberRows,camberCol] = size(camberLine);
    if camberRows ~= 2
        error('Missing coordinates. Cannot interpolate camber line')
    end
    if camberLine(1,end) > 1
        error('Normalize camber line coordinates')
    end
    camberD = zeros(1,length(normControlPoints));
    camber = spline(camberLine(1,:),camberLine(2,:),normControlPoints);
    %Prep Central Diff 
    centDiffCont = zeros(1,2*length(normControlPoints));
    count = 1;
    for i = 1:2:2*length(normControlPoints)
        centDiffCont(i) = normControlPoints(count)-0.01;
        centDiffCont(i+1) = normControlPoints(count)+0.01;
        count = count +1;
    end
    camberCentCont = spline(camberLine(1,:),camberLine(2,:),centDiffCont); 
    %Take central diff
    count = 1;
    for i = 1:2:2*length(normControlPoints)
            camberD(count) = (camberCentCont(i+1)-camberCentCont(i))/0.02;
            count = count + 1;
    end   
    camberDrep = repmat(camberD',2*bPanels,1);
    camberMat  = repmat(camberD',1,2*bPanels);
end

end

