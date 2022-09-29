%VLM 8.0
%Safa Bakhshi
%Set flag to nondimentionalize lift and drag
clear all
NoDimLift = 1;
global timeGlobal
global MAC
global spanGlobal
%Set Planforms 
%   data = [taper;rootC;AR;LambLE(degrees)] column for each section
% Please un-comment the desired planform to run the code
%Simple Rectangular
% sections = 1;
% taper = .1;
% S = 200;
% ARfix = 3;
% b = sqrt(ARfix*S); %Span
% rootC = 3*S/b; %Root Chord
% phi = 0*(pi/180);%Dihedral Angle
% alphaIncRoot = 0*(pi/180); %Incidence at Root converted to radians
% alphaIncTip = 0*(pi/180); %Incidence at Tip converted to radians
% leLambda = 70;
% data = [taper;rootC;ARfix;leLambda];
% camberLine = 0;
% camberLine = '2412';
% camberLine = [0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65...
%  0.7,0.75,0.8,0.85,0.9,0.95,1.;0.000483,0.01105636,0.01848429,0.02444774,0.0287301,0.03171505...
%  0.03373632,0.03501125,0.03544156,0.03533584,0.03463478,0.03314842...
%  0.03116487,0.0287059,0.02578617,0.02241971,0.01859947,0.0144237...
%  0.0099137,0.00509027,0.];

%Simple Delta
% sections =1;
% taper = 0; %Taper Ratio
% S = 100; % Reference Area
% ARfix = 2; %Aspect Ratio
% b = sqrt(ARfix*S); %Span
% rootC = S*2/b; %Root Chord
% alphaIncRoot = 0*(pi/180); %Incidence at Root converted to radians
% alphaIncTip = 0*(pi/180); %Incidence at Tip converted to radians
% leLambda= (90-atand((b/2)/rootC)); %Sweep angle 
% phi = 0*(pi/180); %Dihedral Angle
% data = [taper;rootC;ARfix;leLambda];
% camberLine = 0;

%Simple Delta(Clipped Tips)
sections =1;
taper = 0.031243; %Taper Ratio
S = 98.3882; % Reference Area
ARfix = 1.9077; %Aspect Ratio
b = 13.7002; %Span
rootC = 14.1421; %Root Chord
alphaIncRoot = 0*(pi/180); %Incidence at Root converted to radians
alphaIncTip = 0*(pi/180); %Incidence at Tip converted to radians
leLambda= 63.434948822922000; %Sweep angle 
phi = 0*(pi/180); %Dihedral Angle
data = [taper;rootC;ARfix;leLambda];
camberLine = 0;


%VLM Control: Set panels and modelling options
% bPanels = round(4*b/pi)/2; %Spanwise Panels/2
bPanels = 8;
cPanels = 8; %Chordwise Panels
wakeT = 400; %span lengths for wake termination
% dt = .1*ones(1,100); %Time Step Vector
dt = [.005*ones(1,100) .1*ones(1,100)];
wRelaxT = 1; %Factor of Wake Relaxation
wRelaxS = 1;
% vS = (0.6/2)*dt(1); %Leading edge separation vorticity fraction
vS = 0;
Lfac = 1-vS;

%Freestream Control: Set alpha,beta, and leading edge panels for suction
%analogy, Uinf, and rho
NLE = 0; %Number of rows of panels to apply suction analogy to starting with the Leading Edge
alpha = 20*(pi/180); % Angle of attack in degrees/converted to Radians 
beta = 0*(pi/180); %Sideslip angle in degrees/converted to Radians
Uinf =1; %Freestream velocity 
rho =1; %Freestream Density

%Use Geometry Engine to Obtain Panels
[controlPointsY,controlPointsX,panelQuarterC,panelTQuarterC,panelGeomY,panelGeomX,VcontrolPointsX,S,b,N,K,CAve,geomChord,chord] = geomEngine3(sections,data,bPanels,cPanels,'Plot');

realAR = b^2/S;
spanGlobal = b;
MAC = (1/S)*trapz(panelGeomY,geomChord.^2);
timeGlobal=0.01;
fprintf('The Aspect Ratio of this Wing is :%f\n',realAR);
fprintf('The Mean Aerodynamic Chord is :%f\n',MAC);
fprintf('The Wing Area is: %f\n',S);

%Set Twist
twist = spline([panelGeomY(1) panelGeomY(length(panelGeomY)/2 +0.5) panelGeomY(end)],[alphaIncTip alphaIncRoot alphaIncTip],controlPointsY);
twist = repmat(twist,cPanels,1);
twist = twist(:);
%Set Wing Camber
[camberDrep,camberMat,camber] = camberFit(camberLine,cPanels,bPanels,N);

%VLM Main Process
%Wake Simulation
vorticityMatrix = zeros(N,N);%K1
vorticityWake = zeros(N,N);%Wake Contributions
vorticitySep = zeros(N,N);%Leading Edge Separation Contributions

%Initialize wake
xw = panelGeomX(end,:);
yw = panelGeomY;
zw = zeros(1,length(panelGeomY));
%Initialize wake history 
xwHis = zeros(length(dt)+1,length(xw));
ywHis = zeros(length(dt)+1,length(yw));
zwHis = zeros(length(dt)+1,length(zw));
xwHis(1,:) = xw;
ywHis(1,:) = yw;
zwHis(1,:) = zw;
%Initilize wake velocity field
uw = zeros(1,length(panelGeomY));
vw = zeros(1,length(panelGeomY));
ww = zeros(1,length(panelGeomY));
%Initialize leading edge separation field
xL = panelGeomX(1,:);
yL = panelGeomY;
zL = zeros(1,length(panelGeomY));
%Initialize leading edge separation history 
xLHis = zeros(length(dt)+1,length(xw));
yLHis = zeros(length(dt)+1,length(yw));
zLHis = zeros(length(dt)+1,length(zw));
xLHis(1,:) = xL;
yLHis(1,:) = yL;
zLHis(1,:) = zL;
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
                    if m == 1
                        Tfac = vS;
                    else
                        Tfac = 1;
                    end                   
                    f = timeStep;

                    [u1,v1,w1] = VORTEX(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),0,-1*(panelQuarterC(m,k)-controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                    
                    %Accumulate Trailing Vortex Sheet Contributions
                    %Wake vorticies
                    %Vortex 2: TE to wake
                    %Vortex 3: wake to TE
                    if f ~= 1
                        [u2,v2,w2] = VORTEX(-1*(xwHis(f,k+1) - xwHis(f-1,k+1)),ywHis(f,k+1)-ywHis(f-1,k+1),zwHis(f,k+1)-zwHis(f-1,k+1),...
                            -1*(xwHis(f-1,k+1)-controlPointsX(j,i)),ywHis(f-1,k+1)-controlPointsY(i),zwHis(f-1,k+1));

                        [u3,v3,w3] = VORTEX(-1*(xwHis(f-1,k)- xwHis(f,k)),ywHis(f-1,k)-ywHis(f,k),zwHis(f-1,k)-zwHis(f,k),...
                            -1*(xwHis(f,k)-controlPointsX(j,i)),ywHis(f,k)-controlPointsY(i),zwHis(f,k)); 
                        %Apply effects of camber and dihedral to contributions
                        %before assigning to vorticity matrix. Note
                        %asymmetircal dihedral angle on each wing
                        %Contribution accumulated over timesteps

                       if panelCount > N/2
                            vorticityWake(panelCount,contCount) = vorticityWake(panelCount,contCount) + Tfac*(w2+w3)*cos(phi)*cos(atan(camberMat(j,i))) - Tfac*(u2+u3)*sin(atan(camberMat(j,i)))*cos(phi) - Tfac*(v2+v3)*sin(phi);
                       else
                            vorticityWake(panelCount,contCount) = vorticityWake(panelCount,contCount) + Tfac*(w2+w3)*cos(phi)*cos(atan(camberMat(j,i))) - Tfac*(u2+u3)*sin(atan(camberMat(j,i)))*cos(phi) - Tfac*(v2+v3)*sin(-phi);
                       end
                    end
                    
                    %Accumulate Leading Edge Separated Vorticies
                    %Contributions
                    %Wake vorticies
                    %Vortex 6: TE to wake
                    %Vortex 7: wake to TE
                    if f ~= 1 && m ==1
                        [u6,v6,w6] = VORTEX(-1*(xLHis(f,k+1) - xLHis(f-1,k+1)),yLHis(f,k+1)-yLHis(f-1,k+1),zLHis(f,k+1)-zLHis(f-1,k+1),...
                            -1*(xLHis(f-1,k+1)-controlPointsX(j,i)),yLHis(f-1,k+1)-controlPointsY(i),zLHis(f-1,k+1));

                        [u7,v7,w7] = VORTEX(-1*(xLHis(f-1,k)- xLHis(f,k)),yLHis(f-1,k)-yLHis(f,k),zLHis(f-1,k)-zLHis(f,k),...
                            -1*(xLHis(f,k)-controlPointsX(j,i)),yLHis(f,k)-controlPointsY(i),zLHis(f,k)); 
                        

                        %Apply effects of camber and dihedral to contributions
                        %before assigning to vorticity matrix. Note
                        %asymmetircal dihedral angle on each wing
                        %Contribution accumulated over timesteps
                        
                    
                       if panelCount > N/2
                            vorticitySep(panelCount,contCount) = vorticitySep(panelCount,contCount) + Lfac*(w6+w7)*cos(phi)*cos(atan(camberMat(j,i))) - Lfac*(u6+u7)*sin(atan(camberMat(j,i)))*cos(phi) - Lfac*(v6+v7)*sin(phi);
                       else
                            vorticitySep(panelCount,contCount) = vorticitySep(panelCount,contCount) + Lfac*(w6+w7)*cos(phi)*cos(atan(camberMat(j,i))) - Lfac*(u6+u7)*sin(atan(camberMat(j,i)))*cos(phi) - Lfac*(v6+v7)*sin(-phi);
                       end
                    end                    
                    
                    %Infinite Addon from last known wake points                    
                    
                    [u2inf,v2inf,w2inf] = VORTEX(wakeT*b - (-1*xwHis(f,k+1)),tan(beta)*(wakeT*b - (-1*xwHis(f,k+1))),tan(alpha)*(wakeT*b - (-1*xwHis(f,k+1))),...
                        -1*(xwHis(f,k+1)-controlPointsX(j,i)),ywHis(f,k+1)-controlPointsY(i),zwHis(f,k+1));

                    [u3inf,v3inf,w3inf] = VORTEX((-1*xwHis(f,k))- wakeT*b ,tan(beta)*((-1*xwHis(f,k))- wakeT*b),tan(alpha)*((-1*xwHis(f,k))- wakeT*b),...
                        wakeT*b-(-1*controlPointsX(j,i)),ywHis(f,k)+tan(beta)*(wakeT*b - (-1*xwHis(f,k)))-controlPointsY(i),zwHis(f,k)+tan(alpha)*(wakeT*b-(-1*xwHis(f,k))));
                    if m==1
                    [u6inf,v6inf,w6inf] = VORTEX(wakeT*b - (-1*xLHis(f,k+1)),tan(beta)*(wakeT*b - (-1*xLHis(f,k+1))),tan(alpha)*(wakeT*b - (-1*xLHis(f,k+1))),...
                        -1*(xLHis(f,k+1)-controlPointsX(j,i)),yLHis(f,k+1)-controlPointsY(i),zLHis(f,k+1));

                    [u7inf,v7inf,w7inf] = VORTEX((-1*xLHis(f,k))- wakeT*b ,tan(beta)*((-1*xLHis(f,k))- wakeT*b),tan(alpha)*((-1*xLHis(f,k))- wakeT*b),...
                        wakeT*b-(-1*controlPointsX(j,i)),yLHis(f,k)+tan(beta)*(wakeT*b - (-1*xLHis(f,k)))-controlPointsY(i),zLHis(f,k)+tan(alpha)*(wakeT*b-(-1*xLHis(f,k))));
                    else
                        u6inf =0;v6inf=0;w6inf=0;u7inf =0;v7inf=0;w7inf=0;
                    end
                    u2inf=Tfac*u2inf; v2inf=Tfac*v2inf; w2inf=Tfac*w2inf; u3inf =Tfac*u3inf; v3inf=Tfac*v3inf; w3inf=Tfac*w3inf;
                    u6inf=Lfac*u6inf; v6inf=Lfac*v6inf; w6inf=Lfac*w6inf; u7inf =Lfac*u7inf; v7inf=Lfac*v7inf; w7inf=Lfac*w7inf;
                    %Trailing Vorticies on the wing
                    %Vortex 4: quarter chord to TE
                    %Vortex 5: TE to quarter chord
                    if m == 1
                        [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-controlPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i),0);
                        [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k))),0,0,-1*(panelGeomX(end,k)-controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                        
                        [u4L,v4L,w4L] = VORTEX((-1*(panelGeomX(1,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-controlPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i),0);
                        [u5L,v5L,w5L] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(1,k))),0,0,-1*(panelGeomX(1,k)-controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);                        

                        u4 =Tfac*u4 + Lfac*u4L; v4 =Tfac*v4 + Lfac*v4L; w4=Tfac*w4+ Lfac*w4L; u5=Tfac*u5+ Lfac*u5L; v5=Tfac*v5+ Lfac*v5L; w5=Tfac*w5+ Lfac*w5L;
                    else
                        [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-controlPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i),0);
                        [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k))),0,0,-1*(panelGeomX(end,k)-controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                    end
                    %Sum contributions from each segment in the u,v,w
                    %directions. Excluding the tracked wake.
                    u = u1 + u2inf + u3inf + u6inf + u7inf + u4 + u5;
                    v = v1 + v2inf + v3inf + v6inf + v7inf + v4 + v5; 
                    w = w1 + w2inf + w3inf + w6inf + w7inf + w4 + w5;
                    
                    %Apply effects of camber and dihedral to contributions
                    %before assigning to vorticity matrix. Note
                    %asymmetircal dihedral angle on each wing
                   if panelCount > N/2
                    vorticityMatrix(panelCount,contCount) = vorticityWake(panelCount,contCount) + vorticitySep(panelCount,contCount) + w*cos(phi)*cos(atan(camberMat(j,i))) - u*sin(atan(camberMat(j,i)))*cos(phi) - v*sin(phi);
                   else
                    vorticityMatrix(panelCount,contCount) = vorticityWake(panelCount,contCount) + vorticitySep(panelCount,contCount) + w*cos(phi)*cos(atan(camberMat(j,i))) - u*sin(atan(camberMat(j,i)))*cos(phi) - v*sin(-phi);
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
uL = zeros(1,length(panelGeomY));
vL = zeros(1,length(panelGeomY));
wL = zeros(1,length(panelGeomY));
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
                    u6tot = 0;
                    u7tot = 0;
                    v2tot = 0;
                    v3tot = 0;
                    v6tot = 0;
                    v7tot = 0;
                    w2tot = 0;
                    w3tot = 0;
                    w6tot = 0;
                    w7tot = 0;

                   if m == 1
                        Tfac = vS;
                   else
                        Tfac = 1;
                   end                    
                    for f = 1:timeStep
                        if f ~= 1
                            %Wake Vorticies
                            %Vortex 2: TE to wake
                            %Vortex 3: wake to TE 
                            [u2,v2,w2] = VORTEX(-1*(xwHis(f,k+1) - xwHis(f-1,k+1)),ywHis(f,k+1)-ywHis(f-1,k+1),zwHis(f,k+1)-zwHis(f-1,k+1),...
                                -1*(xwHis(f-1,k+1)-xw(i)),ywHis(f-1,k+1)- yw(i),zwHis(f-1,k+1)-zw(i));

                            [u3,v3,w3] = VORTEX(-1*(xwHis(f-1,k)- xwHis(f,k)),ywHis(f-1,k)-ywHis(f,k),zwHis(f-1,k)-zwHis(f,k),...
                                -1*(xwHis(f,k)-xw(i)),ywHis(f,k)- yw(i),zwHis(f,k)-zw(i));

                            u2 =u2*Tfac; v2 =v2*Tfac; w2=w2*Tfac; u3=u3*Tfac; v3=v3*Tfac; w3=w3*Tfac;
                        
                            if m == 1
                                [u6,v6,w6] = VORTEX(-1*(xLHis(f,k+1) - xLHis(f-1,k+1)),yLHis(f,k+1)-yLHis(f-1,k+1),zLHis(f,k+1)-zLHis(f-1,k+1),...
                                    -1*(xLHis(f-1,k+1)-xw(i)),yLHis(f-1,k+1)- yw(i),zLHis(f-1,k+1)-zw(i));

                                [u7,v7,w7] = VORTEX(-1*(xLHis(f-1,k)- xLHis(f,k)),yLHis(f-1,k)-yLHis(f,k),zLHis(f-1,k)-zLHis(f,k),...
                                    -1*(xLHis(f,k)-xw(i)),yLHis(f,k)- yw(i),zLHis(f,k)-zw(i));  

                                u6 =u6*Lfac; v6 =v6*Lfac; w6 =w6*Lfac; u7=u7*Lfac; v7=v7*Lfac; w7=w7*Lfac;
                            else
                                u6 =0; v6 =0; w6 =0; u7=0; v7=0; w7=0;
                            end
                            u2tot = u2tot + u2;
                            u3tot = u3tot + u3;
                            u6tot = u6tot + u6;
                            u7tot = u7tot + u7;  
                            v2tot = v2tot + v2;
                            v3tot = v3tot + v3;
                            v6tot = v6tot + v6;
                            v7tot = v7tot + v7;
                            w2tot = w2tot + w2;
                            w3tot = w3tot + w3;
                            w6tot = w6tot + w6;
                            w7tot = w7tot + w7;

                        end
                        if f == timeStep
                            %Infinite Addon from last known wake points
                            [u2inf,v2inf,w2inf] = VORTEX(wakeT*b - (-1*xwHis(f,k+1)),tan(beta)*(wakeT*b - (-1*xwHis(f,k+1))),tan(alpha)*(wakeT*b - (-1*xwHis(f,k+1))),...
                                -1*(xwHis(f,k+1)-xw(i)),ywHis(f,k+1)-yw(i),zwHis(f,k+1)-zw(i));

                            [u3inf,v3inf,w3inf] = VORTEX((-1*xwHis(f,k))- wakeT*b ,tan(beta)*((-1*xwHis(f,k))- wakeT*b),tan(alpha)*((-1*xwHis(f,k))- wakeT*b),...
                                wakeT*b-(-1*xw(i)),ywHis(f,k)+tan(beta)*(wakeT*b-(-1*xwHis(f,k)))-yw(i),zwHis(f,k)+tan(alpha)*(wakeT*b-(-1*xwHis(f,k)))-zw(i));
                            if m==1
                                [u6inf,v6inf,w6inf] = VORTEX(wakeT*b - (-1*xLHis(f,k+1)),tan(beta)*(wakeT*b - (-1*xLHis(f,k+1))),tan(alpha)*(wakeT*b - (-1*xLHis(f,k+1))),...
                                    -1*(xLHis(f,k+1)-xw(i)),yLHis(f,k+1)-yw(i),zLHis(f,k+1)-zw(i));

                                [u7inf,v7inf,w7inf] = VORTEX((-1*xLHis(f,k))- wakeT*b ,tan(beta)*((-1*xLHis(f,k))- wakeT*b),tan(alpha)*((-1*xLHis(f,k))- wakeT*b),...
                                    wakeT*b-(-1*xw(i)),yLHis(f,k)+tan(beta)*(wakeT*b-(-1*xLHis(f,k)))-yw(i),zLHis(f,k)+tan(alpha)*(wakeT*b-(-1*xLHis(f,k)))-zw(i));
                            else
                                u6inf =0;v6inf=0;w6inf=0;u7inf =0;v7inf=0;w7inf=0;
                            end
                            u2inf=Tfac*u2inf; v2inf=Tfac*v2inf; w2inf=Tfac*w2inf; u3inf =Tfac*u3inf; v3inf=Tfac*v3inf; w3inf=Tfac*w3inf;  
                            u6inf=Lfac*u6inf; v6inf=Lfac*v6inf; w6inf=Lfac*w6inf; u7inf =Lfac*u7inf; v7inf=Lfac*v7inf; w7inf=Lfac*w7inf;  
                            
                            u2tot = u2tot + u2inf;
                            u3tot = u3tot + u3inf;
                            u6tot = u6tot + u6inf;
                            u7tot = u7tot + u7inf;
                            v2tot = v2tot + v2inf;
                            v3tot = v3tot + v3inf;
                            v6tot = v6tot + v6inf;
                            v7tot = v7tot + v7inf;
                            w2tot = w2tot + w2inf;
                            w3tot = w3tot + w3inf;
                            w6tot = w6tot + w6inf;
                            w7tot = w7tot + w7inf;
                        end

                    end
                    %Trailing Vorticies on the wing
                    %Vortex 4: quarter chord to TE
                    %Vortex 5: TE to quarter chord
                    if m == 1
                        [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-xw(i)),panelGeomY(k+1)-yw(i),-zw(i));
                        [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k))),0,0,-1*(panelGeomX(end,k)-xw(i)),panelGeomY(k)-yw(i),-zw(i));
                        
                        [u4L,v4L,w4L] = VORTEX((-1*(panelGeomX(1,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-xw(i)),panelGeomY(k+1)-yw(i),-zw(i));
                        [u5L,v5L,w5L] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(1,k))),0,0,-1*(panelGeomX(1,k)-xw(i)),panelGeomY(k)-yw(i),-zw(i));                        

                        u4 =Tfac*u4 + Lfac*u4L; v4 =Tfac*v4 + Lfac*v4L; w4=Tfac*w4+ Lfac*w4L; u5=Tfac*u5+ Lfac*u5L; v5=Tfac*v5+ Lfac*v5L; w5=Tfac*w5+ Lfac*w5L;
                    else
                        [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-xw(i)),panelGeomY(k+1)-yw(i),-zw(i));
                        [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k))),0,0,-1*(panelGeomX(end,k)-xw(i)),panelGeomY(k)-yw(i),-zw(i));
                    end
                    %Sum contributions from each segment in the u,v,w
                    %directions. Include factor of relaxation here

                    uw(i) = uw(i) + wRelaxS*gammaMatrix(m,k)*(u1 + u2tot + u3tot + u4 + u5 + u6tot + u7tot);
                    vw(i) = vw(i) + wRelaxS*gammaMatrix(m,k)*(v1 + v2tot + v3tot + v4 + v5 + v6tot + v7tot); 
                    ww(i) = ww(i) + wRelaxS*gammaMatrix(m,k)*(w1 + w2tot + w3tot + w4 + w5 + w6tot + w7tot);      

        end  
    end 
end

%Loop through each tracked leading edge vortex
for i = 1:length(panelGeomY)
    %Loop through each panel on the wing. Calculate downwash contributions BY this
    %panel 
        for k = 1:length(controlPointsY)
            for m = 1:cPanels
                %Bound Vortex: Vortex 1

                [u1,v1,w1] = VORTEX(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),0,-1*(panelQuarterC(m,k)-xL(i)),panelGeomY(k)-yL(i),-zL(i));

                %Trailing Vortex Contributions. Recalculated on each
                %time step
                u2tot = 0;
                u3tot = 0;
                u6tot = 0;
                u7tot = 0;
                v2tot = 0;
                v3tot = 0;
                v6tot = 0;
                v7tot = 0;
                w2tot = 0;
                w3tot = 0;
                w6tot = 0;
                w7tot = 0;
                
               if m == 1
                    Tfac = vS;
               else
                    Tfac = 1;
               end          
                for f = 1:timeStep
                    if f ~= 1
                        %Wake Vorticies
                        %Vortex 2: TE to wake
                        %Vortex 3: wake to TE                            

                        [u2,v2,w2] = VORTEX(-1*(xwHis(f,k+1) - xwHis(f-1,k+1)),ywHis(f,k+1)-ywHis(f-1,k+1),zwHis(f,k+1)-zwHis(f-1,k+1),...
                            -1*(xwHis(f-1,k+1)-xL(i)),ywHis(f-1,k+1)- yL(i),zwHis(f-1,k+1)-zL(i));

                        [u3,v3,w3] = VORTEX(-1*(xwHis(f-1,k)- xwHis(f,k)),ywHis(f-1,k)-ywHis(f,k),zwHis(f-1,k)-zwHis(f,k),...
                            -1*(xwHis(f,k)-xL(i)),ywHis(f,k)- yL(i),zwHis(f,k)-zL(i));       

                        u2 =u2*Tfac; v2 =v2*Tfac; w2=w2*Tfac; u3=u3*Tfac; v3=v3*Tfac; w3=w3*Tfac;


                        if m==1
                            [u6,v6,w6] = VORTEX(-1*(xLHis(f,k+1) - xLHis(f-1,k+1)),yLHis(f,k+1)-yLHis(f-1,k+1),zLHis(f,k+1)-zLHis(f-1,k+1),...
                                -1*(xLHis(f-1,k+1)-xL(i)),yLHis(f-1,k+1)- yL(i),zLHis(f-1,k+1)-zL(i));

                            [u7,v7,w7] = VORTEX(-1*(xLHis(f-1,k)- xLHis(f,k)),yLHis(f-1,k)-yLHis(f,k),zLHis(f-1,k)-zLHis(f,k),...
                                -1*(xLHis(f,k)-xL(i)),yLHis(f,k)- yL(i),zLHis(f,k)-zL(i));      

                            u6 =u6*Lfac; v6 =v6*Lfac; w6 =w6*Lfac; u7=u7*Lfac; v7=v7*Lfac; w7=w7*Lfac;
                        end
                            u2tot = u2tot + u2;
                            u3tot = u3tot + u3;
                            u6tot = u6tot + u6;
                            u7tot = u7tot + u7;
                            v2tot = v2tot + v2;
                            v3tot = v3tot + v3;
                            v6tot = v6tot + v6;
                            v7tot = v7tot + v7;
                            w2tot = w2tot + w2;
                            w3tot = w3tot + w3;
                            w6tot = w6tot + w6;
                            w7tot = w7tot + w7;
                    end
                    if f == timeStep
                        %Infinite Addon from last known wake points
                        [u2inf,v2inf,w2inf] = VORTEX(wakeT*b - (-1*xwHis(f,k+1)),tan(beta)*(wakeT*b - (-1*xwHis(f,k+1))),tan(alpha)*(wakeT*b - (-1*xwHis(f,k+1))),...
                            -1*(xwHis(f,k+1)-xL(i)),ywHis(f,k+1)-yL(i),zwHis(f,k+1)-zL(i));

                        [u3inf,v3inf,w3inf] = VORTEX((-1*xwHis(f,k))- wakeT*b ,tan(beta)*((-1*xwHis(f,k))- wakeT*b),tan(alpha)*((-1*xwHis(f,k))- wakeT*b),...
                            wakeT*b-(-1*xL(i)),ywHis(f,k)+tan(beta)*(wakeT*b-(-1*xwHis(f,k)))-yL(i),zwHis(f,k)+tan(alpha)*(wakeT*b-(-1*xwHis(f,k)))-zL(i));
                        if m ==1
                            [u6inf,v6inf,w6inf] = VORTEX(wakeT*b - (-1*xLHis(f,k+1)),tan(beta)*(wakeT*b - (-1*xLHis(f,k+1))),tan(alpha)*(wakeT*b - (-1*xLHis(f,k+1))),...
                                -1*(xLHis(f,k+1)-xL(i)),yLHis(f,k+1)-yL(i),zLHis(f,k+1)-zL(i));

                            [u7inf,v7inf,w7inf] = VORTEX((-1*xLHis(f,k))- wakeT*b ,tan(beta)*((-1*xLHis(f,k))- wakeT*b),tan(alpha)*((-1*xLHis(f,k))- wakeT*b),...
                                wakeT*b-(-1*xL(i)),yLHis(f,k)+tan(beta)*(wakeT*b-(-1*xLHis(f,k)))-yL(i),zLHis(f,k)+tan(alpha)*(wakeT*b-(-1*xLHis(f,k)))-zL(i));
                        else
                            u6inf =0;v6inf=0;w6inf=0;u7inf =0;v7inf=0;w7inf=0;
                        end
                        u2inf=Tfac*u2inf; v2inf=Tfac*v2inf; w2inf=Tfac*w2inf; u3inf =Tfac*u3inf; v3inf=Tfac*v3inf; w3inf=Tfac*w3inf;  
                        u6inf=Lfac*u6inf; v6inf=Lfac*v6inf; w6inf=Lfac*w6inf; u7inf =Lfac*u7inf; v7inf=Lfac*v7inf; w7inf=Lfac*w7inf;
                        
                        u2tot = u2tot + u2inf;
                        u3tot = u3tot + u3inf;
                        u6tot = u6tot + u6inf;
                        u7tot = u7tot + u7inf;
                        v2tot = v2tot + v2inf;
                        v3tot = v3tot + v3inf;
                        v6tot = v6tot + v6inf;
                        v7tot = v7tot + v7inf;
                        w2tot = w2tot + w2inf;
                        w3tot = w3tot + w3inf;
                        w6tot = w6tot + w6inf;
                        w7tot = w7tot + w7inf;
                    end

                end
                %Trailing Vorticies on the wing
                %Vortex 4: quarter chord to TE
                %Vortex 5: TE to quarter chord
                if m == 1
                    [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-xL(i)),panelGeomY(k+1)-yL(i),-zL(i));
                    [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k))),0,0,-1*(panelGeomX(end,k)-xL(i)),panelGeomY(k)-yL(i),-zL(i));

                    [u4L,v4L,w4L] = VORTEX((-1*(panelGeomX(1,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-xL(i)),panelGeomY(k+1)-yL(i),-zL(i));
                    [u5L,v5L,w5L] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(1,k))),0,0,-1*(panelGeomX(1,k)-xL(i)),panelGeomY(k)-yL(i),-zL(i));                        

                    u4 =Tfac*u4 + Lfac*u4L; v4 =Tfac*v4 + Lfac*v4L; w4=Tfac*w4+ Lfac*w4L; u5=Tfac*u5+ Lfac*u5L; v5=Tfac*v5+ Lfac*v5L; w5=Tfac*w5+ Lfac*w5L;
                else
                    [u4,v4,w4] = VORTEX((-1*(panelGeomX(end,k+1)) - (-1*panelQuarterC(m,k+1))),0,0,-1*(panelQuarterC(m,k+1)-xL(i)),panelGeomY(k+1)-yL(i),-zL(i));
                    [u5,v5,w5] = VORTEX((-1*panelQuarterC(m,k))- (-1*(panelGeomX(end,k))),0,0,-1*(panelGeomX(end,k)-xL(i)),panelGeomY(k)-yL(i),-zL(i));
                end
                %Sum contributions from each segment in the u,v,w
                %directions. Include factor of relaxation here
                uL(i) = uL(i) + wRelaxT*gammaMatrix(m,k)*(u1 + u2tot + u3tot + u4 + u5 + u6tot + u7tot);
                vL(i) = vL(i) + wRelaxT*gammaMatrix(m,k)*(v1 + v2tot + v3tot + v4 + v5 + v6tot + v7tot); 
                wL(i) = wL(i) + wRelaxT*gammaMatrix(m,k)*(w1 + w2tot + w3tot + w4 + w5 + w6tot + w7tot);

            end           
        end
end  

%Use Forward Euler to advance one time step and obtain the K1 term
xwK1 = -dt(timeStep)*(uw + Uinf*cos(alpha)*cos(beta));
ywK1 = dt(timeStep)*(vw + Uinf*cos(alpha)*sin(beta));
zwK1 = dt(timeStep)*(ww + Uinf*sin(alpha));
xLK1 = -dt(timeStep)*(uL + Uinf*cos(alpha)*cos(beta));
yLK1 = dt(timeStep)*(vL + Uinf*cos(alpha)*sin(beta));
zLK1 = dt(timeStep)*(wL + Uinf*sin(alpha));

%Solver 1: Forward Euler
% Apply solver to track next wake location
xw = xw + (xwK1);
yw = yw + (ywK1);
zw = zw + (zwK1);
xL = xL + xLK1;
yL = yL + yLK1;
zL = zL + zLK1;
%Solver 2: RK2
% Apply solver to track next wake location
% xw = xw + (xwK1+xwK2)/2;
% yw = yw + (ywK1+ywK2)/2;
% zw = zw + (zwK1+zwK2)/2;
%Update wake position history
xwHis(timeStep+1,:) = xw;
ywHis(timeStep+1,:) = yw;
zwHis(timeStep+1,:) = zw;
xLHis(timeStep+1,:) = xL;
yLHis(timeStep+1,:) = yL;
zLHis(timeStep+1,:) = zL;
timeGlobal = timeGlobal + dt(timeStep);
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

%Plotting Section

%Plot Wing Descretization with horseshoe vortices and wake
wakeTPlot = wakeT*b/100;
drawSim3(controlPointsY,controlPointsX,panelQuarterC,panelGeomY,panelGeomX,alpha,beta,wakeTPlot,xwHis,ywHis,zwHis,xLHis,yLHis,zLHis,dt,bPanels,cPanels,MAC,'freeWake','NOanimate');

%Plot Panel Geometry and Control Points
%Create a Y coordinate matrix for the contour plot
bigY = zeros(cPanels,2*sum(bPanels));
for n = 1:cPanels
    bigY(n,:) = controlPointsY;
end

%Plot a total lift contour plot
figure()
contourf(bigY,VcontrolPointsX,potentialLiftDist,90,'linecolor','none')
colorbar
title('Total lift distribution of planform');
xlabel('y');
ylabel('-x');
colormap jet

%Calculate spanwise lift and drag distributions by summing each column of the
%lift and drag contours
potliftDistribution = sum(potentialLiftDist);


%Divide by K to obtain the lift per unit span
for t = 1:sections
          spanwiseTotalLiftDist([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])  = ...
              potliftDistribution([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])/K(t);
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


%Calculate and Print Coefficients
LiftPot = sum(potliftDistribution);
CLPot = 2*LiftPot/(rho*S*Uinf^2);
fprintf('CL = %f\n',CLPot);

CDiTreffz = 2*inducedDragTreffz/(rho*Uinf^2*S);
fprintf('CDi Treffz %f\n',CDiTreffz);

einv =  CLPot^2/(pi*realAR)/CDiTreffz;
fprintf('Pot e = %f\n',einv);



% function [u,v,w] = VORTEX0(Gx, Gy, Gz, Rx, Ry, Rz)
%     % Compute velocity induced by unit strength vortex filament.
%     % x is positive downstream, y is positive from centerline to right tip, z is upward.
%     % Gx, Gy, Gz are the lengths of the vortex filament in the x, y, and z directions.
%     % Rx, Ry, Rz is a vector from the vortex root to the control point: Rx = xroot-xctl
%     R = [Rx,Ry,Rz]; G = [Gx,Gy,Gz]; R2 = dot(R,R); G2 = dot(G,G); eps = 1E-2;
%     GR = dot(G,R); GXR = cross(G,R); GXR2 = dot(GXR,GXR); GPR = G+R; GPR2 = dot(GPR,GPR);
%     Vtot = (GR/sqrt(R2)-(G2+GR)/sqrt(GPR2))/(4.0*pi*GXR2+eps); V = Vtot*GXR;
%     u=V(1);v=V(2);w=V(3);
%   end


%Function to calculate 
function [u,v,w] = VORTEX(Gx, Gy, Gz, Rx, Ry, Rz)
%Vortex Root at Origin
x1n = 0;
y1n = 0;
z1n = 0;

%Vortex Tip Location
x2n = Gx;
y2n = Gy;
z2n = Gz;

% Set Control Point location

%Convention 2: Kroo convention(Vortex Root minus Control Point)
x = -Rx;
y = -Ry;
z = -Rz;

% 
G2=Gx*Gx+Gy*Gy+Gz*Gz;
R2=Rx*Rx+Ry*Ry+Rz*Rz;
E2=(Gx+Rx)*(Gx+Rx)+(Gy+Ry)*(Gy+Ry)+(Gz+Rz)*(Gz+Rz);
R1Length = sqrt(R2);
GLength = sqrt(G2);
R2Length = sqrt(E2);
spp = (R1Length+GLength+R2Length)/2;
h = 4*spp*(spp-R1Length)*(spp-GLength)*(spp-R2Length)/G2;

% h = norm(cross([Rx,Ry,Rz],[Gx,Gy,Gz]))/norm([Gx,Gy,Gz]);
MAC = 14.1421;
if h < (1E-3*MAC)^2
% if h < sqrt(4*.001*timeGlobal)
    u=0.0;
    v=0.0;
    w=0.0;
else
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

    V = (sig1-sig2) / (4.0*pi);
    u = V*(R1XR2x/R1XR2);
    v = V*(R1XR2y/R1XR2);
    w = V*(R1XR2z/R1XR2);
end



end