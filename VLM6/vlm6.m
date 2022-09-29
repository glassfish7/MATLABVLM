%VLM 6.0
%Safa Bakhshi
%Set flag to nondimentionalize lift and drag
NoDimLift = 1;

%Set Planforms 
%   data = [taper;rootC;AR;LambLE(degrees)] column for each section
% Please un-comment the desired planform to run the code

%Simple Rectangular
% sections = 1;
% taper = 1;
% S = 100;
% ARfix = 10;
% b = sqrt(ARfix*S); %Span
% rootC = S/b; %Root Chord
% phi = 0*(pi/180);%Dihedral Angle
% alphaIncRoot = 0*(pi/180); %Incidence at Root converted to radians
% alphaIncTip = 0*(pi/180); %Incidence at Tip converted to radians
% leLambda = 0;
% data = [taper;rootC;ARfix;leLambda];
% camberLine = 0;
% camberLine = '2412';
% camberLine = [0.0468750000000000,0.109375000000000,0.171875000000000,0.234375000000000,0.296875000000000,0.359375000000000,0.421875000000000,...
% 0.484375000000000,0.546875000000000,0.609375000000000,0.671875000000000,0.734375000000000,0.796875000000000,0.859375000000000,0.921875000000000,...
% 0.984375000000000;0.0130723741319444,0.0153076171875000,0.0171088324652778,0.0184760199652778,0.0194091796875000,0.0199083116319444,0.0199734157986111,...
% 0.0196044921875000,0.0188015407986111,0.0175645616319444,0.0158935546875000,0.0137885199652778,0.0112494574652778,0.00827636718750000,0.00486924913194445,...
% 0.00102810329861111];

% sections = 1;
% taper = 0.8181;
% S = 30;
% ARfix = 7.5;
% b = sqrt(ARfix*S); %Span
% rootC = 2.2; %Root Chord
% phi = 0*(pi/180);%Dihedral Angle
% alphaIncRoot = 2*(pi/180); %Incidence at Root converted to radians
% alphaIncTip = 2*(pi/180); %Incidence at Tip converted to radians
% leLambda = 0;
% data = [taper;rootC;ARfix;leLambda];
% camberLine = 0;
%Simple Delta
sections =1;
taper = 0; %Taper Ratio
S = 100; % Reference Area
ARfix = 2; %Aspect Ratio
b = sqrt(ARfix*S); %Span
rootC = S*2/b; %Root Chord
alphaIncRoot = 0*(pi/180); %Incidence at Root converted to radians
alphaIncTip = 0*(pi/180); %Incidence at Tip converted to radians
leLambda= (90-atand((b/2)/rootC)); %Sweep angle 
phi = 0*(pi/180); %Dihedral Angle
data = [taper;rootC;ARfix;leLambda];
camberLine = 0;

%VLM Control: Set panels and modelling options
% bPanels = round(4*b/pi)/2; %Spanwise Panels/2
bPanels= 8;
cPanels = 8; %Chordwise Panels
wakeT = 400; %span lengths for wake termination
os = 0; %Kink Trailing Edge Offset Factor
dt = 0.1*ones(1,1); %Time Step for wake progresion
wRelax = 1; %Factor of Wake Relaxation

%Freestream Control: Set alpha,beta, and leading edge panels for suction
%analogy, Uinf, and rho
NLE = 0; %Number of rows of panels to apply suction analogy to starting with the Leading Edge
alpha = 20*(pi/180); % Angle of attack in degrees/converted to Radians 
beta = 0*(pi/180); %Sideslip angle in degrees/converted to Radians
Uinf =1; %Freestream velocity 
rho =1; %Freestream Density

%Use Geometry Engine to Obtain Panels
[controlPointsY,controlPointsX,panelQuarterC,panelGeomY,panelGeomX,VcontrolPointsX,S,b,N,K,CAve,~,chord] = geomEngine2(sections,data,bPanels,cPanels,'NoPlot');

realAR = b^2/S;
fprintf('The Aspect Ratio of this Wing is :%f\n',realAR);
fprintf('The Wing Area is: %f\n',S);

%Set Twist
twist = spline([panelGeomY(1) panelGeomY(length(panelGeomY)/2 +0.5) panelGeomY(end)],[alphaIncTip alphaIncRoot alphaIncTip],controlPointsY);
twist = repmat(twist,cPanels,1);
twist = twist(:);
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

%VLM Main Process
%Wake Simulation
vorticityMatrix = zeros(N,N);%K1
vorticityWake = zeros(N,N);%K1
vorticityMatrixK2 = zeros(N,N);%K2
vorticityWakeK2 = zeros(N,N);%K2

%Initialize wake
xw = panelGeomX(end,:);
yw = panelGeomY;
zw = zeros(1,length(panelGeomY));
%Initialize wake history 
xwHis = zeros(length(dt),length(xw));
ywHis = zeros(length(dt),length(yw));
zwHis = zeros(length(dt),length(zw));
xwHis(1,:) = xw;
ywHis(1,:) = yw;
zwHis(1,:) = zw;
%Initilize wake velocity field
uw = zeros(1,length(panelGeomY));
vw = zeros(1,length(panelGeomY));
ww = zeros(1,length(panelGeomY));
%Initilize progress bar
bar = waitbar(0,'Time Progress');

%Time Step Process
for timeStep = 1:length(dt)
waitbar(timeStep/length(dt),bar)%Update progress bar

%Predictor Step: K1
%Set counters to 1
panelCount =1;
contCount =1;
%Loop through each panel on the wing. Calculate downwash contributions ONTO this
%panel 
for i = 1:length(controlPointsY)
    for j = 1:cPanels
        %Loop through each panel on the wing. Calculate downwash contributions BY this
        %panel 
        for k = 1:length(controlPointsY)  
           for m = 1:cPanels
                    f = timeStep;
                    %Accumulate Trailing Vortex Contributions
                    if f ~= 1
                        %Wake vorticies
                        %Vortex 2: TE to wake
                        %Vortex 3: wake to TE
                        [u2,v2,w2] = VORTEX(-1*(xwHis(f,k+1) - xwHis(f-1,k+1)),ywHis(f,k+1)-ywHis(f-1,k+1),zwHis(f,k+1)-zwHis(f-1,k+1),...
                            -1*(xwHis(f-1,k+1)-controlPointsX(j,i)),ywHis(f-1,k+1)-controlPointsY(i),zwHis(f-1,k+1));

                        [u3,v3,w3] = VORTEX(-1*(xwHis(f-1,k)- xwHis(f,k)),ywHis(f-1,k)-ywHis(f,k),zwHis(f-1,k)-zwHis(f,k),...
                            -1*(xwHis(f,k)-controlPointsX(j,i)),ywHis(f,k)-controlPointsY(i),zwHis(f,k)); 
 
                        %Apply effects of camber and dihedral to contributions
                        %before assigning to vorticity matrix. Note
                        %asymmetircal dihedral angle on each wing
                        %Contribution accumulated over timesteps
                       if panelCount > N/2
                            vorticityWake(panelCount,contCount) = vorticityWake(panelCount,contCount) + (w2+w3)*cos(phi)*cos(atan(camberMat(j,i))) - (u2+u3)*sin(atan(camberMat(j,i)))*cos(phi) - (v2+v3)*sin(phi);
                       else
                            vorticityWake(panelCount,contCount) = vorticityWake(panelCount,contCount) + (w2+w3)*cos(phi)*cos(atan(camberMat(j,i))) - (u2+u3)*sin(atan(camberMat(j,i)))*cos(phi) - (v2+v3)*sin(-phi);
                       end
                    end
                    %Infinite Addon from last known wake points
                    [u2inf,v2inf,w2inf] = VORTEX(wakeT*b - (-1*xwHis(f,k+1)),tan(beta)*(wakeT*b - (-1*xwHis(f,k+1))),tan(alpha)*(wakeT*b - (-1*xwHis(f,k+1))),...
                        -1*(xwHis(f,k+1)-controlPointsX(j,i)),ywHis(f,k+1)-controlPointsY(i),zwHis(f,k+1));

                    [u3inf,v3inf,w3inf] = VORTEX((-1*xwHis(f,k))- wakeT*b ,tan(beta)*((-1*xwHis(f,k))- wakeT*b),tan(alpha)*((-1*xwHis(f,k))- wakeT*b),...
                        wakeT*b-(-1*controlPointsX(j,i)),ywHis(f,k)+tan(beta)*(wakeT*b - (-1*xwHis(f,k)))-controlPointsY(i),zwHis(f,k)+tan(alpha)*(wakeT*b-(-1*xwHis(f,k))));
                    %Trailing Vorticies on the wing
                    %Vortex 4: quarter chord to TE
                    %Vortex 5: TE to quarter chord
                    [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-controlPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i),0);
                    [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k))),0,0,-1*(panelGeomX(end,k)-controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                    

                    %Bound Vortex: Vortex 1
                    [u1,v1,w1] = VORTEX(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),0,-1*(panelQuarterC(m,k)-controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);


                    %Sum contributions from each segment in the u,v,w
                    %directions. Excluding the tracked wake.
                    u = u1 + u2inf + u3inf + u4 + u5;
                    v = v1 + v2inf + v3inf + v4 + v5; 
                    w = w1 + w2inf + w3inf + w4 + w5;
                    
                    %Apply effects of camber and dihedral to contributions
                    %before assigning to vorticity matrix. Note
                    %asymmetircal dihedral angle on each wing
                   if panelCount > N/2
                    vorticityMatrix(panelCount,contCount) = vorticityWake(panelCount,contCount) + w*cos(phi)*cos(atan(camberMat(j,i))) - u*sin(atan(camberMat(j,i)))*cos(phi) - v*sin(phi);
                   else
                    vorticityMatrix(panelCount,contCount) = vorticityWake(panelCount,contCount) + w*cos(phi)*cos(atan(camberMat(j,i))) - u*sin(atan(camberMat(j,i)))*cos(phi) - v*sin(-phi);
                   end
                   contCount = contCount + 1;    
            end
        end      
        %Update counters
        panelCount = panelCount + 1;
        contCount =1;
    end
end

%Generate RHS of matrix equation with no penetration boundary condition.
%Includes effects from alpha, beta, phi, and camber
%RHS1 is left wing and RHS2 is right wing
%Concatinate vectors 
RHS1 = Uinf*(sin(atan(camberDrep(1:N/2,:)))*cos(phi)*cos(alpha)*cos(beta)+...
    sin(beta)*sin(-phi)...
    -sin(alpha)*cos(beta)*cos(atan(camberDrep(1:N/2,:)))*cos(phi)).*ones(N/2,1);
RHS2 = Uinf*(sin(atan(camberDrep(N/2+1:end,:)))*cos(phi)*cos(alpha)*cos(beta)+...
    sin(beta)*sin(phi)-...
    sin(alpha)*cos(beta)*cos(atan(camberDrep(N/2+1:end,:)))*cos(phi)).*ones(N/2,1);
RHS = [RHS1;RHS2];
%No dihedral version for debug
% RHS = (Uinf*(-sin(alpha+twist)+camberDrep)).*ones(N,1);

%Solve for vorticity
%Alternative Approach: (1/sqrt(n) * R) \ ((1/sqrt(n)  * Q' * (dA * V)))
gamma = linsolve(vorticityMatrix,RHS);

%Place gamma vector back into array
gammaMatrix = zeros(cPanels,2*sum(bPanels));
gammaCount =1;
for p = 1:2*sum(bPanels)
    for o = 1:cPanels
        gammaMatrix(o,p) = gamma(gammaCount);
        gammaCount = gammaCount + 1;
    end
end
%Wake Tracking Process
%Reset wake velocity field purturbations
uw = zeros(1,length(panelGeomY));
vw = zeros(1,length(panelGeomY));
ww = zeros(1,length(panelGeomY));
%Loop through each tracked line vortex
for i = 1:length(panelGeomY)
        %Loop through each panel on the wing. Calculate downwash contributions BY this
        %panel 
        for k = 1:length(controlPointsY)
            for m = 1:cPanels
                    %Bound Vortex: Vortex 1
                    [u1,v1,w1] = VORTEX(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),0,-1*(panelQuarterC(m,k)-xw(i)),panelGeomY(k)-yw(i),-zw(i));
                    %Trailing Vortex Contributions. Recalculated on each
                    %time step
                    u2tot = 0;
                    u3tot = 0;
                    v2tot = 0;
                    v3tot = 0;
                    w2tot = 0;
                    w3tot = 0;
                    for f = 1:timeStep
                        if f ~= 1
                            %Wake Vorticies
                            %Vortex 2: TE to wake
                            %Vortex 3: wake to TE
                            [u2,v2,w2] = VORTEX(-1*(xwHis(f,k+1) - xwHis(f-1,k+1)),ywHis(f,k+1)-ywHis(f-1,k+1),zwHis(f,k+1)-zwHis(f-1,k+1),...
                                -1*(xwHis(f-1,k+1)-xw(i)),ywHis(f-1,k+1)- yw(i),zwHis(f-1,k+1)-zw(i));

                            [u3,v3,w3] = VORTEX(-1*(xwHis(f-1,k)- xwHis(f,k)),ywHis(f-1,k)-ywHis(f,k),zwHis(f-1,k)-zwHis(f,k),...
                                -1*(xwHis(f,k)-xw(i)),ywHis(f,k)- yw(i),zwHis(f,k)-zw(i));
                            u2tot = u2tot + u2;
                            u3tot = u3tot + u3;
                            v2tot = v2tot + v2;
                            v3tot = v3tot + v3;
                            w2tot = w2tot + w2;
                            w3tot = w3tot + w3;
                        end
                        if f == timeStep
                            %Infinite Addon from last known wake points
                            [u2inf,v2inf,w2inf] = VORTEX(wakeT*b - (-1*xwHis(f,k+1)),tan(beta)*(wakeT*b - (-1*xwHis(f,k+1))),tan(alpha)*(wakeT*b - (-1*xwHis(f,k+1))),...
                                -1*(xwHis(f,k+1)-xw(i)),ywHis(f,k+1)-yw(i),zwHis(f,k+1)-zw(i));

                            [u3inf,v3inf,w3inf] = VORTEX((-1*xwHis(f,k))- wakeT*b ,tan(beta)*((-1*xwHis(f,k))- wakeT*b),tan(alpha)*((-1*xwHis(f,k))- wakeT*b),...
                                wakeT*b-(-1*xw(i)),ywHis(f,k)+tan(beta)*(wakeT*b-(-1*xwHis(f,k)))-yw(i),zwHis(f,k)+tan(alpha)*(wakeT*b-(-1*xwHis(f,k)))-zw(i));
                            u2tot = u2tot + u2inf;
                            u3tot = u3tot + u3inf;
                            v2tot = v2tot + v2inf;
                            v3tot = v3tot + v3inf;
                            w2tot = w2tot + w2inf;
                            w3tot = w3tot + w3inf;                       
                        end
                    end
                    %Trailing Vorticies on the wing
                    %Vortex 4: quarter chord to TE
                    %Vortex 5: TE to quarter chord
                    [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-xw(i)),panelGeomY(k+1)-yw(i),-zw(i));
                    [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k))),0,0,-1*(panelGeomX(end,k)-xw(i)),panelGeomY(k)-yw(i),-zw(i));
                    
                    %Sum contributions from each segment in the u,v,w
                    %directions. Include factor of relaxation here
                    uw(i) = uw(i) + wRelax*gammaMatrix(m,k)*(u1 + u2tot + u3tot + u4 + u5);
                    vw(i) = vw(i) + wRelax*gammaMatrix(m,k)*(v1 + v2tot + v3tot + v4 + v5); 
                    ww(i) = ww(i) + wRelax*gammaMatrix(m,k)*(w1 + w2tot + w3tot + w4 + w5);
            end
        end
end
%Use Forward Euler to advance one time step and obtain the K1 term
xwK1 = -dt(timeStep)*(uw + Uinf*cos(alpha)*cos(beta));
ywK1 = dt(timeStep)*(vw + Uinf*cos(alpha)*sin(beta));
zwK1 = dt(timeStep)*(ww + Uinf*sin(alpha));

%Corrector Step: K2
%Set counters to 1
panelCount =1;
contCount =1;
%Loop through each panel on the wing. Calculate downwash contributions ONTO this
%panel 
for i = 1:length(controlPointsY)
    for j = 1:cPanels
        %Loop through each panel on the wing. Calculate downwash contributions BY this
        %panel 
        for k = 1:length(controlPointsY)
            for m = 1:cPanels
                f = timeStep;
                    %Trailing Vortex Contributions. Accumulated over each
                    %time step
                    if f ~= 1
                        %Wake vorticies
                        %Vortex 2: TE to wake
                        %Vortex 3: wake to TE
                        [u2,v2,w2] = VORTEX(-1*((xwHis(f,k+1)+xwK1(k+1)) - xwHis(f-1,k+1)),ywHis(f,k+1)+ywK1(k+1)-ywHis(f-1,k+1),zwHis(f,k+1)+zwK1(k+1)-zwHis(f-1,k+1),...
                            -1*(xwHis(f-1,k+1)-controlPointsX(j,i)),ywHis(f-1,k+1)-controlPointsY(i),zwHis(f-1,k+1));

                        [u3,v3,w3] = VORTEX(-1*(xwHis(f-1,k)- (xwHis(f,k)+xwK1(k))),ywHis(f-1,k)-(ywHis(f,k)+ywK1(k)),zwHis(f-1,k)-(zwHis(f,k)+zwK1(k)),...
                            -1*((xwHis(f,k)+xwK1(k))-controlPointsX(j,i)),(ywHis(f,k)+ywK1(k))-controlPointsY(i),(zwHis(f,k)+zwK1(k)));         
                        %Apply effects of camber and dihedral to contributions
                        %before assigning to vorticity matrix. Note
                        %asymmetircal dihedral angle on each wing
                        %Contribution accumulated over each time step
                       if panelCount > N/2
                            vorticityWakeK2(panelCount,contCount) = vorticityWakeK2(panelCount,contCount) + (w2+w3)*cos(phi)*cos(atan(camberMat(j,i))) - (u2+u3)*sin(atan(camberMat(j,i)))*cos(phi) - (v2+v3)*sin(phi);
                       else
                            vorticityWakeK2(panelCount,contCount) = vorticityWakeK2(panelCount,contCount) + (w2+w3)*cos(phi)*cos(atan(camberMat(j,i))) - (u2+u3)*sin(atan(camberMat(j,i)))*cos(phi) - (v2+v3)*sin(-phi);
                       end
                    end
                    %Infinite Addon from last known wake points
                    [u2inf,v2inf,w2inf] = VORTEX(wakeT*b - (-1*(xwHis(f,k+1)+xwK1(k+1))),tan(beta)*(wakeT*b - (-1*(xwHis(f,k+1)+xwK1(k+1)))),tan(alpha)*(wakeT*b - (-1*(xwHis(f,k+1)+xwK1(k+1)))),...
                        -1*((xwHis(f,k+1)+xwK1(k+1))-controlPointsX(j,i)),ywHis(f,k+1)+ywK1(k+1)-controlPointsY(i),zwHis(f,k+1)+zwK1(k+1));

                    [u3inf,v3inf,w3inf] = VORTEX((-1*(xwHis(f,k)+xwK1(k)))- wakeT*b ,tan(beta)*((-1*(xwHis(f,k)+xwK1(k)))- wakeT*b),tan(alpha)*((-1*(xwHis(f,k)+xwK1(k)))- wakeT*b),...
                        wakeT*b-(-1*controlPointsX(j,i)),(ywHis(f,k)+ywK1(k))+tan(beta)*(wakeT*b - (-1*(xwHis(f,k)+xwK1(k))))-controlPointsY(i),(zwHis(f,k)+zwK1(k))+tan(alpha)*(wakeT*b-(-1*(xwHis(f,k)+xwK1(k)))));
                    %Trailing Vorticies on the wing
                    %Vortex 4: quarter chord to TE
                    %Vortex 5: TE to quarter chord
                    [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-controlPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i),0);
                    [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k))),0,0,-1*(panelGeomX(end,k)-controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                    %Bound Vortex: Vortex 1
                    [u1,v1,w1] = VORTEX(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),0,-1*(panelQuarterC(m,k)-controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                    %Sum contributions from each segment in the u,v,w
                    %directions
                    u = u1 + u2inf + u3inf + u4 + u5;
                    v = v1 + v2inf + v3inf + v4 + v5; 
                    w = w1 + w2inf + w3inf + w4 + w5;
                    %Apply effects of camber and dihedral to contributions
                    %before assigning to vorticity matrix. Note
                    %asymmetircal dihedral angle on each wing
                   if panelCount > N/2
                    vorticityMatrixK2(panelCount,contCount) = vorticityWakeK2(panelCount,contCount) + w*cos(phi)*cos(atan(camberMat(j,i))) - u*sin(atan(camberMat(j,i)))*cos(phi) - v*sin(phi);
                   else
                    vorticityMatrixK2(panelCount,contCount) = vorticityWakeK2(panelCount,contCount) + w*cos(phi)*cos(atan(camberMat(j,i))) - u*sin(atan(camberMat(j,i)))*cos(phi) - v*sin(-phi);
                   end
                   contCount = contCount + 1;    
            end
        end      
        %Update counters
        panelCount = panelCount + 1;
        contCount =1;
    end
end

%Generate RHS of matrix equation with no penetration boundary condition.
%Includes effects from alpha, beta, phi, and camber
%RHS1 is left wing and RHS2 is right wing
%Concatinate vectors 
RHS1 = Uinf*(sin(atan(camberDrep(1:N/2,:)))*cos(phi)*cos(alpha)*cos(beta)+...
    sin(beta)*sin(-phi)...
    -sin(alpha)*cos(beta)*cos(atan(camberDrep(1:N/2,:)))*cos(phi)).*ones(N/2,1);
RHS2 = Uinf*(sin(atan(camberDrep(N/2+1:end,:)))*cos(phi)*cos(alpha)*cos(beta)+...
    sin(beta)*sin(phi)-...
    sin(alpha)*cos(beta)*cos(atan(camberDrep(N/2+1:end,:)))*cos(phi)).*ones(N/2,1);
RHS = [RHS1;RHS2];
%No dihedral version for debug
% RHS = (Uinf*(-sin(alpha+twist)+camberDrep)).*ones(N,1);

%Solve for vorticity
gammaK2 = linsolve(vorticityMatrixK2,RHS);

%Place gamma vector back into array
gammaMatrixK2 = zeros(cPanels,2*sum(bPanels));
gammaCount =1;
for p = 1:2*sum(bPanels)
    for o = 1:cPanels
        gammaMatrixK2(o,p) = gammaK2(gammaCount);
        gammaCount = gammaCount + 1;
    end
end

%Reset wake velocity field purturbations
uw = zeros(1,length(panelGeomY));
vw = zeros(1,length(panelGeomY));
ww = zeros(1,length(panelGeomY));
%Repeat wake velocity field calculations for f(x+K1) in order to obtain K2
%for the RK2 solver
for i = 1:length(panelGeomY)
        %Loop through each panel on the wing. Calculate downwash contributions BY this
        %panel 
        for k = 1:length(controlPointsY)
            for m = 1:cPanels
                    %Bound Vortex: Vortex 1
                    [u1,v1,w1] = VORTEX(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),0,-1*(panelQuarterC(m,k)-(xw(i)+xwK1(i))),panelGeomY(k)-(yw(i)+ywK1(i)),-(zw(i)+zwK1(i)));
                    %Trailing Vortex Contributions. Recalculated on each
                    %time step
                    u2tot = 0;
                    u3tot = 0;
                    v2tot = 0;
                    v3tot = 0;
                    w2tot = 0;
                    w3tot = 0;
                    for f = 1:timeStep
                        if f ~= 1
                        %Wake Vorticies
                        %Vortex 2: TE to wake
                        %Vortex 3: wake to TE
                        [u2,v2,w2] = VORTEX(-1*((xwHis(f,k+1)+xwK1(k+1)) - xwHis(f-1,k+1)),ywHis(f,k+1)+ywK1(k+1)-ywHis(f-1,k+1),zwHis(f,k+1)+zwK1(k+1)-zwHis(f-1,k+1),...
                            -1*(xwHis(f-1,k+1)-(xw(i)+xwK1(i))),ywHis(f-1,k+1)-(yw(i)+ywK1(i)),zwHis(f-1,k+1)-(zw(i)+zwK1(i)));

                        [u3,v3,w3] = VORTEX(-1*(xwHis(f-1,k)- (xwHis(f,k)+xwK1(k))),ywHis(f-1,k)-(ywHis(f,k)+ywK1(k)),zwHis(f-1,k)-(zwHis(f,k)+zwK1(k)),...
                            -1*((xwHis(f,k)+xwK1(k))-(xw(i)+xwK1(i))),(ywHis(f,k)+ywK1(k))-(yw(i)+ywK1(i)),(zwHis(f,k)+zwK1(k))-(zw(i)+zwK1(i)));  
                        u2tot = u2tot + u2;
                        u3tot = u3tot + u3;
                        v2tot = v2tot + v2;
                        v3tot = v3tot + v3;
                        w2tot = w2tot + w2;
                        w3tot = w3tot + w3;
                        end
                        if f == timeStep
                            %Infinite Addon from last known wake points
                            [u2inf,v2inf,w2inf] = VORTEX(wakeT*b - (-1*(xwHis(f,k+1)+xwK1(k+1))),tan(beta)*(wakeT*b - (-1*(xwHis(f,k+1)+xwK1(k+1)))),tan(alpha)*(wakeT*b - (-1*(xwHis(f,k+1)+xwK1(k+1)))),...
                                -1*((xwHis(f,k+1)+xwK1(k+1))-(xw(i)+xwK1(i))),ywHis(f,k+1)+ywK1(k+1)-(yw(i)+ywK1(i)),zwHis(f,k+1)+zwK1(k+1)-(zw(i)+zwK1(i)));

                            [u3inf,v3inf,w3inf] = VORTEX((-1*(xwHis(f,k)+xwK1(k)))- wakeT*b ,tan(beta)*((-1*(xwHis(f,k)+xwK1(k)))- wakeT*b),tan(alpha)*((-1*(xwHis(f,k)+xwK1(k)))- wakeT*b),...
                                wakeT*b-(-1*(xw(i)+xwK1(i))),(ywHis(f,k)+ywK1(k))+tan(beta)*(wakeT*b - (-1*(xwHis(f,k)+xwK1(k))))-(yw(i)+ywK1(i)),(zwHis(f,k)+zwK1(k))+tan(alpha)*(wakeT*b-(-1*(xwHis(f,k)+xwK1(k))))-(zw(i)+zwK1(i)));
                            u2tot = u2tot + u2inf;
                            u3tot = u3tot + u3inf;
                            v2tot = v2tot + v2inf;
                            v3tot = v3tot + v3inf;
                            w2tot = w2tot + w2inf;
                            w3tot = w3tot + w3inf;                       
                        end
                    end
                    %Trailing Vorticies on the wing
                    %Vortex 4: quarter chord to TE
                    %Vortex 5: TE to quarter chord
                    [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-(xw(i)+xwK1(i))),panelGeomY(k+1)-(yw(i)+ywK1(i)),-(zw(i)+zwK1(i)));
                    [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k))),0,0,-1*(panelGeomX(end,k)-(xw(i)+xwK1(i))),panelGeomY(k)-(yw(i)+ywK1(i)),-(zw(i)+zwK1(i)));
                    
                    %Sum contributions from each segment in the u,v,w
                    %directions. Include factor of relaxation here
                    uw(i) = uw(i) + wRelax*gammaMatrixK2(m,k)*(u1 + u2tot + u3tot + u4 + u5);
                    vw(i) = vw(i) + wRelax*gammaMatrixK2(m,k)*(v1 + v2tot + v3tot + v4 + v5); 
                    ww(i) = ww(i) + wRelax*gammaMatrixK2(m,k)*(w1 + w2tot + w3tot + w4 + w5);
            end
        end
end
%Advance one time step and obtain the K2 term
xwK2 = -dt(timeStep)*(uw + Uinf*cos(alpha)*cos(beta));
ywK2 = dt(timeStep)*(vw + Uinf*cos(alpha)*sin(beta));
zwK2 = dt(timeStep)*(ww + Uinf*sin(alpha));
%Solver 1: Forward Euler
% Apply solver to track next wake location
% xw = xw + (xwK1);
% yw = yw + (ywK1);
% zw = zw + (zwK1);
%Solver 2: RK2
% Apply solver to track next wake location
xw = xw + (xwK1+xwK2)/2;
yw = yw + (ywK1+ywK2)/2;
zw = zw + (zwK1+zwK2)/2;
%Update wake position history
xwHis(timeStep+1,:) = xw;
ywHis(timeStep+1,:) = yw;
zwHis(timeStep+1,:) = zw;
end

%This section of code breaks the gammaMatrix back into sections 
%by identifying start and ends points of each section in the matrix
sPoints = zeros(2,sections);
ePoints = zeros(2,sections);
for i = 1:sections
    sPoints(1,i) = 1 + sum(bPanels(sections:-1:i))-bPanels(i);
    sPoints(2,i) = sPoints(1,i) + sum(bPanels(1:i))+sum(bPanels(1:i-1));

    ePoints(1,i) = sPoints(1,i) + bPanels(i)-1;
    ePoints(2,i) = sPoints(2,i) + bPanels(i)-1;
end
%Use the start and end points to apply Kutta-Joukowski to each section
%with the appropriate width(K) of each section
for t = 1:sections
          potentialLiftDist(:,[sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])  = ...
              rho*Uinf*gammaMatrix(:,[sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])*K(t);

end

% In order to apply the Kutta-Joukowski the bound vortices are projected
% normal to the freestream. NOTE: Sweep is already accounted for in the panel width K from
% geomEngine.m
for cPan = 1:cPanels
    potentialLiftDist(cPan,1:bPanels) = potentialLiftDist(cPan,1:bPanels)*cos(-beta);
    potentialLiftDist(cPan,bPanels+1:end) = potentialLiftDist(cPan,bPanels+1:end)*cos(beta);
end

%Calculate the induced drag produced by the wing in potential flow with Treffz Plane Analysis
gammaTreffz = sum(gammaMatrix);
wiTreffz = zeros(1,length(controlPointsY));
%Find downwash at the wake
for i = 1:length(controlPointsY)
                    wiTreffz(i) = sum(gammaTreffz.*(-1./(panelGeomY(1:end-1)-controlPointsY(i)) + 1./(panelGeomY(2:end)-controlPointsY(i))))/(4*pi);   
end
inducedDragTreffz = (rho)*sum(gammaTreffz.*wiTreffz);
CDiTreffz = 2*inducedDragTreffz/(rho*Uinf^2*S);
fprintf('CDi Treffz %f\n',CDiTreffz);

%Vortex Lift Section
%Calculate Zero Leading edge suction lift 
%Find Normal Force on each panel for zero-suction lift
%Also calculate the potential induced drag
normalForce = zeros(cPanels,2*sum(bPanels));
potIndDrag =  zeros(cPanels,2*sum(bPanels));
suctionNormalForce = zeros(cPanels,2*sum(bPanels));
panelCount = 1;

for i = 1:length(controlPointsY)
    for j = 1:cPanels
        wi = 0;
        %Calculate downwash on panel      
        for k = 1:length(controlPointsY)
            for m = 1:cPanels         
                    %Bound Vortex: Vortex 1
                    [u1,v1,w1] = VORTEX(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),0,-1*(panelQuarterC(m,k)-VcontrolPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                    
                    %Trailing Vortex Contributions
                    u2tot = 0;
                    u3tot = 0;
                    v2tot = 0;
                    v3tot = 0;
                    w2tot = 0;
                    w3tot = 0;
                    for f = 1:length(dt)
                        if f ~= 1
                            %Wake vorticies
                            %Vortex 2: TE to wake
                            %Vortex 3: wake to TE
                            [u2,v2,w2] = VORTEX(-1*(xwHis(f,k+1) - xwHis(f-1,k+1)),ywHis(f,k+1)-ywHis(f-1,k+1),zwHis(f,k+1)-zwHis(f-1,k+1),...
                                -1*(xwHis(f-1,k+1)-VcontrolPointsX(j,i)),ywHis(f-1,k+1)-controlPointsY(i),zwHis(f-1,k+1));

                            [u3,v3,w3] = VORTEX(-1*(xwHis(f-1,k)- xwHis(f,k)),ywHis(f-1,k)-ywHis(f,k),zwHis(f-1,k)-zwHis(f,k),...
                                -1*(xwHis(f,k)-VcontrolPointsX(j,i)),ywHis(f,k)-controlPointsY(i),zwHis(f,k));                           
                            u2tot = u2tot + u2;
                            u3tot = u3tot + u3;
                            v2tot = v2tot + v2;
                            v3tot = v3tot + v3;
                            w2tot = w2tot + w2;
                            w3tot = w3tot + w3; 
                        end
                        
                        if f == length(dt)
                            %Infinite Addon from last known wake points
                            [u2inf,v2inf,w2inf] = VORTEX(wakeT*b - (-1*xwHis(f,k+1)),tan(beta)*(wakeT*b - (-1*xwHis(f,k+1))),tan(alpha)*(wakeT*b - (-1*xwHis(f,k+1))),...
                                -1*(xwHis(f,k+1)-VcontrolPointsX(j,i)),ywHis(f,k+1)-controlPointsY(i),zwHis(f,k+1));

                            [u3inf,v3inf,w3inf] = VORTEX((-1*xwHis(f,k))- wakeT*b ,tan(beta)*((-1*xwHis(f,k))- wakeT*b),tan(alpha)*((-1*xwHis(f,k))- wakeT*b),...
                                wakeT*b-(-1*VcontrolPointsX(j,i)),ywHis(f,k)+tan(beta)*(wakeT*b - (-1*xwHis(f,k)))-controlPointsY(i),zwHis(f,k)+tan(alpha)*(wakeT*b-(-1*xwHis(f,k))));
                            u2tot = u2tot + u2inf;
                            u3tot = u3tot + u3inf;
                            v2tot = v2tot + v2inf;
                            v3tot = v3tot + v3inf;
                            w2tot = w2tot + w2inf;
                            w3tot = w3tot + w3inf; 
                        end
                    end
                    %Trailing Vorticies on the wing
                    %Vortex 4: quarter chord to TE
                    %Vortex 5: TE to quarter chord
                    [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-VcontrolPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i),0);
                    [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k))),0,0,-1*(panelGeomX(end,k)-VcontrolPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                    
                    %Sum contributions from each segment in the w
                    %direction
                    w = w1 + w2tot + w3tot + w4 + w5;
                    
                    %Apply effects of camber and dihedral to contributions
                    %before assigning to vorticity matrix. Note
                    %asymmetircal dihedral angle on each wing
                   if panelCount > N/2
                    wi = wi + gammaMatrix(m,k)*w*cos(phi)*cos(atan(camberMat(j,i)));
                   else
                    wi = wi + gammaMatrix(m,k)*w*cos(phi)*cos(atan(camberMat(j,i)));
                   end        
            end
        end 
       alphai = atan(-wi/Uinf);
       inducedDragi = potentialLiftDist(j,i)*sin(alphai);
       potIndDrag(j,i) = inducedDragi;
       if j <= NLE
           %Calculate Lift due to suction if panels part of LE
           normalForce(j,i) = potentialLiftDist(j,i)*cos(alpha) + inducedDragi*sin(alpha);
           dxLE = panelGeomX(j,i+1) - panelGeomX(j,i);
           dy = panelGeomY(i+1) - panelGeomY(i);
           angleLE = atan(dxLE/dy);
           suctionNormalForce(j,i) = (normalForce(j,i)*sin(alpha)- inducedDragi)/cos(angleLE);
       end
       panelCount = panelCount + 1;
    end
end
%Calculate Vortex Lift and Drag
vortexLift =(suctionNormalForce+normalForce)*cos(alpha);
vortexInducedDrag = vortexLift*tan(alpha);
%Calculate No Suction Lift and Drag
noSucLift = (normalForce)*cos(alpha);
noSucInducedDrag = normalForce*tan(alpha);
%Potential Lift and Drag on remaining planform
remainingPlanformLift = zeros(cPanels,2*sum(bPanels));
remainingPlanformLift(NLE+1:end,:) =  potentialLiftDist(NLE+1:end,:);
remainingPlanformDrag = zeros(cPanels,2*sum(bPanels));
remainingPlanformDrag(NLE+1:end,:) =  potIndDrag(NLE+1:end,:);
%Combine the remaining planform lift and drag with the leading edge Vortex
%Lift and Drag
combineLift = vortexLift + remainingPlanformLift;
combineDrag = vortexInducedDrag + remainingPlanformDrag;
%Combine the remaining planform lift and drag with the leading edge No
%Suction Lift and Drag
combineDragNV = noSucInducedDrag + remainingPlanformDrag;
combineLiftNV = noSucLift + remainingPlanformLift;

%Plotting Section

%Plot Wing Descretization with horseshoe vortices and wake
wakeTPlot = wakeT*b/100;
drawSim(controlPointsY,controlPointsX,panelQuarterC,panelGeomY,panelGeomX,alpha,beta,wakeTPlot,xwHis,ywHis,zwHis,dt,bPanels,cPanels,'realWake','NOanimate')

%Plot Panel Geometry and Control Points
%Create a Y coordinate matrix for the contour plot
bigY = zeros(cPanels,2*sum(bPanels));
for n = 1:cPanels
    bigY(n,:) = controlPointsY;
end

%Plot a total lift contour plot
figure()
contourf(bigY,VcontrolPointsX,combineLift,90,'linecolor','none')
colorbar
title('Total lift distribution of planform');
xlabel('y');
ylabel('-x');
colormap jet

%Calculate spanwise lift and drag distributions by summing each column of the
%lift and drag contours
totalLiftDistribution = sum(combineLift);
potliftDistribution = sum(potentialLiftDist);
NVLiftDistribution = sum(combineLiftNV);
potIndDragDistribution = sum(potIndDrag);
totalIndDragDistribution = sum(combineDrag);
NVDist = sum(combineDragNV);

%Divide by K to obtain the lift per unit span
for t = 1:sections
          spanwiseTotalLiftDist([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])  = ...
              totalLiftDistribution([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])/K(t);
          
          spanwiseTotalIndDragDistribution([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])  = ...
              totalIndDragDistribution([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])/K(t);
end
% Plot spanwise total Lift distribution
figure()
if NoDimLift == 1
    plot([-b/2 controlPointsY b/2],[0 spanwiseTotalLiftDist/(.5*Uinf^2*rho*CAve) 0])
    title('Calculated spanwise Total Lift distribution');
    xlabel('y')
    ylabel('L''/q C_ave');
else
    plot([-b/2 controlPointsY b/2],[0 spanwiseTotalLiftDist 0])
    title('Calculated spanwise Total Lift distribution');
    xlabel('y')
    ylabel('L''');
end

% Plot spanwise total induced Drag distribution
figure()
if NoDimLift == 1
    plot([-b/2 controlPointsY b/2],[0 spanwiseTotalIndDragDistribution/(.5*Uinf^2*rho*CAve) 0])
    title('Calculated spanwise total induced drag distribution');
    xlabel('y')
    ylabel('D''/q C_ave');
else
    plot([-b/2 controlPointsY b/2],[0 spanwiseTotalIndDragDistribution 0])
    title('Calculated spanwise total induced drag distribution');
    xlabel('y')
    ylabel('D''');
end

%Calculate and Print Coefficients
LiftPot = sum(potliftDistribution);
CLPot = 2*LiftPot/(rho*S*Uinf^2);
fprintf('CL Potential Lift = %f\n',CLPot);

TotalLift = sum(totalLiftDistribution);
CL = 2*TotalLift/(rho*S*Uinf^2);
fprintf('CL Vortex Lift = %f\n',CL);

NoSucLift = sum(NVLiftDistribution);
CLNV = 2*NoSucLift/(rho*S*Uinf^2);
fprintf('CL No Suction = %f\n',CLNV);

PotentialIndDrag = sum(potIndDragDistribution);
CDiPot = 2*PotentialIndDrag/(rho*S*Uinf^2);
fprintf('CDi = %f\n',CDiPot);

IndDrag = sum(totalIndDragDistribution);
CDiTotal = 2*IndDrag/(rho*S*Uinf^2);
fprintf('CDi Vortex = %f\n',CDiTotal);

NVIndDrag = sum(NVDist);
CDiNV = 2*NVIndDrag/(rho*S*Uinf^2);
fprintf('CDi No Suction = %f\n',CDiNV);

einv =  CLPot^2/(pi*realAR)/CDiPot;
fprintf('Pot e = %f\n',einv);

%Function to calculate 
function [u,v,w] = VORTEX(Gx, Gy, Gz, Rx, Ry, Rz)
    % Compute velocity induced by unit strength vortex filament.
    % x is positive downstream, y is positive from centerline to right tip, z is upward.
    % Gx, Gy, Gz are the lengths of the vortex filament in the x, y, and z directions.
    % Rx, Ry, Rz is a vector from the vortex root to the control point: Rx = xroot-xctl
    R = [Rx,Ry,Rz]; G = [Gx,Gy,Gz]; R2 = dot(R,R); G2 = dot(G,G); eps = 1E-8;
    GR = dot(G,R); GXR = cross(G,R); GXR2 = dot(GXR,GXR); GPR = G+R; GPR2 = dot(GPR,GPR);
    Vtot = (GR/sqrt(R2)-(G2+GR)/sqrt(GPR2))/(4.0*pi*GXR2+eps); V = Vtot*GXR;
    u=V(1);v=V(2);w=V(3);
end
function [u,v,w] = VORTEX0(Gx, Gy, Gz, Rx, Ry, Rz)

%Vortex Root at Origin
x1n = 0;
y1n = 0;
z1n = 0;

%Vortex Tip Location
x2n = Gx;
y2n = Gy;
z2n = Gz;

%Set Control Point location

%Convention 1: Textbook convention(Control Point minus Vortex Root)
% x = Rx;
% y = Ry;
% z = Rz;

%Convention 2: Kroo convention(Vortex Root minus Control Point)
x = -Rx;
y = -Ry;
z = -Rz;
%R1 X R2
R1XR2x = 1*((y-y1n)*(z-z2n)-(y-y2n)*(z-z1n));
R1XR2y = -1*((x-x1n)*(z-z2n)-(x-x2n)*(z-z1n));
R1XR2z = 1*((x-x1n)*(y-y2n)-(x-x2n)*(y-y1n));

%|R1 X R2|
R1XR2 = R1XR2x^2 + R1XR2y^2 + R1XR2z^2;

%sigma: see Virginia Tech textbook ch 6 pg. 18
sig1 = ( (x2n-x1n)*(x-x1n) + (y2n-y1n)*(y-y1n) + (z2n-z1n)*(z-z1n))/...
    sqrt((x-x1n)^2 + (y-y1n)^2 + (z-z1n)^2);
sig2 = ( (x2n-x1n)*(x-x2n) + (y2n-y1n)*(y-y2n) + (z2n-z1n)*(z-z2n))/...
    sqrt((x-x2n)^2 + (y-y2n)^2 + (z-z2n)^2);

%Check if vortex is colinear with the control point
TOL = 1.0E-10;
G2=Gx*Gx+Gy*Gy+Gz*Gz;
R2=Rx*Rx+Ry*Ry+Rz*Rz;
TOL2 = TOL * G2;
 
GXRx=Gy*Rz-Gz*Ry;
GXRy=Gz*Rx-Gx*Rz;
GXRz=Gx*Ry-Gy*Rx;

E1=GXRx*GXRx+GXRy*GXRy+GXRz*GXRz;
E2=(Gx+Rx)*(Gx+Rx)+(Gy+Ry)*(Gy+Ry)+(Gz+Rz)*(Gz+Rz);
ERROR = 0;
if R2 <= TOL2 
    ERROR = 1; 
end
if E2 <= TOL2 
    ERROR = 1; 
end
if E1 <= TOL2*R2 
    ERROR = 1; 
end

if ERROR==0
    V = (sig1-sig2) / (4.0*pi);
    u = V*(R1XR2x/R1XR2);
    v = V*(R1XR2y/R1XR2);
    w = V*(R1XR2z/R1XR2);
else
    u=0.0;
    v=0.0;
    w=0.0;
end
end