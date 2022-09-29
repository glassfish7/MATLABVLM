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
bPanels = 1;
cPanels = 2; %Chordwise Panels
wakeT = 400; %span lengths for wake termination
dt = .01*ones(1,100); %Time Step Vector
% dt = [.005*ones(1,100) .1*ones(1,100)];
wRelaxT = 1; %Factor of Wake Relaxation
wRelaxS = 1;
% vS = (0.6/2)*dt(1); %Leading edge separation vorticity fraction
vS = 0.5;
Lfac = 1-vS;

%Freestream Control: Set alpha,beta, and leading edge panels for suction
%analogy, Uinf, and rho
NLE = 0; %Number of rows of panels to apply suction analogy to starting with the Leading Edge
alpha = 20*(pi/180); % Angle of attack in degrees/converted to Radians 
beta = 0*(pi/180); %Sideslip angle in degrees/converted to Radians
Uinf =1; %Freestream velocity 
rho =1; %Freestream Density

%Use Geometry Engine to Obtain Panels
[controlPointsY,controlPointsX,panelQuarterC,panelQuarterCX,panelTQuarterC,panelGeomY,panelGeomX,VcontrolPointsX,VcontrolPointsY,S,b,N,K,CAve,geomChord,chord] = geomEngine3Conc(sections,data,bPanels,cPanels,'Plot');

realAR = b^2/S;
spanGlobal = b;
% MAC = (1/S)*trapz(panelGeomY,geomChord.^2);
MAC = 9.4412;
timeGlobal=0.01;
fprintf('The Aspect Ratio of this Wing is :%f\n',realAR);
fprintf('The Mean Aerodynamic Chord is :%f\n',MAC);
fprintf('The Wing Area is: %f\n',S);

%Set Twist
% twist = spline([panelGeomY(1) panelGeomY(length(panelGeomY)/2 +0.5) panelGeomY(end)],[alphaIncTip alphaIncRoot alphaIncTip],controlPointsY);
twist=0;
twist = repmat(twist,cPanels,1);
twist = twist(:);
%Set Wing Camber
[camberDrep,camberMat,camber] = camberFit(camberLine,cPanels,bPanels,N);

%VLM Main Process
%Wake Simulation
vorticityMatrix = zeros(N,N);%K1
vorticityWake = zeros(N,N);%K1
vorticitySep = zeros(N,N);%Leading Edge Separation

%Initialize wake
xwInit = zeros(1,2*bPanels+1);
ywInit = panelGeomY(end,:);
zwInit = zeros(1,2*bPanels+1);
xwHis = zeros(1,2*bPanels+1);
ywHis = panelGeomY(end,:);
zwHis = zeros(1,2*bPanels+1);
%Initialize wake history 
% xwHis = zeros(length(dt)+1,length(xw));
% ywHis = zeros(length(dt)+1,length(yw));
% zwHis = zeros(length(dt)+1,length(zw));
% xwHis(1,:) = xw;
% ywHis(1,:) = yw;
% zwHis(1,:) = zw;
%Initilize wake velocity field
uw = zeros(1,2*bPanels+1);
vw = zeros(1,length(panelGeomY(end,:)));
ww = zeros(1,2*bPanels+1);
%Initialize leading edge separation field
xLInit = [panelQuarterCX(end:-1:1,1)' panelQuarterCX(:,1)'];
yLInit = [panelQuarterC(end:-1:1,1)' panelQuarterC(:,end)'];
zLInit = zeros(1,length(xLInit));
xLHis = [panelQuarterCX(end:-1:1,1)' panelQuarterCX(:,1)'];
yLHis = [panelQuarterC(end:-1:1,1)' panelQuarterC(:,end)'];
zLHis = zeros(1,length(xLHis));
%Leading edge velocity field
uL = zeros(1,length([panelQuarterCX(end:-1:1,1)' panelQuarterCX(:,1)']));
vL = zeros(1,length([panelQuarterC(end:-1:1,1)' panelQuarterC(:,end)']));
wL = zeros(1,length(uL));

%Initialize leading edge separation history 
% xLHis = zeros(length(dt)+1,length(xL));
% yLHis = zeros(length(dt)+1,length(yL));
% zLHis = zeros(length(dt)+1,length(zL));
% xLHis(1,:) = xL;
% yLHis(1,:) = yL;
% zLHis(1,:) = zL;
%Initilize progress bar
bar = waitbar(0,'Time Progress');

%Time Step Process
for timeStep = 1:length(dt)
waitbar(timeStep/length(dt),bar)%Update progress bar
vorticityWake = zeros(N,N);%K1
vorticitySep = zeros(N,N);%Leading Edge Separation


uw = zeros(timeStep,2*bPanels+1);
vw = zeros(timeStep,length(panelGeomY(end,:)));
ww = zeros(timeStep,2*bPanels+1);
uL = zeros(timeStep,length([panelQuarterCX(end:-1:1,1)' panelQuarterCX(:,1)']));
vL = zeros(timeStep,length([panelQuarterC(end:-1:1,1)' panelQuarterC(:,end)']));
wL = zeros(timeStep,length([panelQuarterCX(end:-1:1,1)' panelQuarterCX(:,1)']));


%Predictor Step: K1
%Set counters to 1
panelCount =1;
contCount =1;
%Loop through each panel on the wing. Calculate downwash contributions ONTO this
%panel 
for i = 1:2*bPanels
    for j = 1:cPanels
        %Loop through each panel on the wing. Calculate downwash contributions BY this
        %panel 
        for k = 1:2*bPanels
                for m = 1:cPanels
                    if k == 1||k==2*bPanels
                        Tfac = vS;
                    else
                        Tfac = 1;
                    end                   
%                     f = timeStep;

                    [u1,v1,w1] = VORTEX(0,panelQuarterC(m,k+1)-panelQuarterC(m,k),0,-1*(panelQuarterCX(m,k)-controlPointsX(j,i)),panelQuarterC(m,k)-controlPointsY(j,i),0);


                    %Accumulate Trailing Vortex Sheet Contributions
                    %Wake vorticies
                    %Vortex 2: TE to wake
                    %Vortex 3: wake to TE
                    for f = 1:timeStep
                    if f ~= 1
                        [u2,v2,w2] = VORTEX(-1*(xwHis(f,k+1) - xwHis(f-1,k+1)),ywHis(f,k+1)-ywHis(f-1,k+1),zwHis(f,k+1)-zwHis(f-1,k+1),...
                            -1*(xwHis(f-1,k+1)-controlPointsX(j,i)),ywHis(f-1,k+1)-controlPointsY(j,i),zwHis(f-1,k+1));

                        [u3,v3,w3] = VORTEX(-1*(xwHis(f-1,k)- xwHis(f,k)),ywHis(f-1,k)-ywHis(f,k),zwHis(f-1,k)-zwHis(f,k),...
                            -1*(xwHis(f,k)-controlPointsX(j,i)),ywHis(f,k)-controlPointsY(j,i),zwHis(f,k)); 
                        
                        if k==2*bPanels
                            u2 = u2*Tfac; v2=v2*Tfac; w2=w2*Tfac;
                        elseif k==1
                            u3 = u3*Tfac; v3=v3*Tfac; w3=w3*Tfac;
                        end
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
                    
                    %Accumulate Leading Edge Separated Vorticies
                    %Contributions
                    %Wake vorticies
                    %Vortex 6: TE to wake
                    %Vortex 7: wake to TE
                    if f ~= 1
                        if k == 2*bPanels
                         LeadSel = m+cPanels;
                        [u6,v6,w6] = VORTEX(-1*(xLHis(f,LeadSel) - xLHis(f-1,LeadSel)),yLHis(f,LeadSel)-yLHis(f-1,LeadSel),zLHis(f,LeadSel)-zLHis(f-1,LeadSel),...
                            -1*(xLHis(f-1,LeadSel)-controlPointsX(j,i)),yLHis(f-1,LeadSel)-controlPointsY(j,i),zLHis(f-1,LeadSel));
                        u7=0;v7=0;w7=0;
                        u6=u6*Lfac;v6=v6*Lfac;w6=w6*Lfac;
                        elseif k == 1
                         LeadSel = m;
                        [u7,v7,w7] = VORTEX(-1*(xLHis(f-1,LeadSel)- xLHis(f,LeadSel)),yLHis(f-1,LeadSel)-yLHis(f,LeadSel),zLHis(f-1,LeadSel)-zLHis(f,LeadSel),...
                            -1*(xLHis(f,LeadSel)-controlPointsX(j,i)),yLHis(f,LeadSel)-controlPointsY(j,i),zLHis(f,LeadSel)); 
                        u6=0;v6=0;w6=0;
                        u7=u7*Lfac;v7=v7*Lfac;w7=w7*Lfac;
                        end
                        %Apply effects of camber and dihedral to contributions
                        %before assigning to vorticity matrix. Note
                        %asymmetircal dihedral angle on each wing
                        %Contribution accumulated over timesteps
                        
                    
                       if panelCount > N/2
                            vorticitySep(panelCount,contCount) = vorticitySep(panelCount,contCount) + (w6+w7)*cos(phi)*cos(atan(camberMat(j,i))) - (u6+u7)*sin(atan(camberMat(j,i)))*cos(phi) - (v6+v7)*sin(phi);
                       else
                            vorticitySep(panelCount,contCount) = vorticitySep(panelCount,contCount) + (w6+w7)*cos(phi)*cos(atan(camberMat(j,i))) - (u6+u7)*sin(atan(camberMat(j,i)))*cos(phi) - (v6+v7)*sin(-phi);
                       end
                    end                    
                    
                    %Infinite Addon from last known wake points                    
                    if f==timeStep
                    [u2inf,v2inf,w2inf] = VORTEX(wakeT*b - (-1*xwHis(f,k+1)),tan(beta)*(wakeT*b - (-1*xwHis(f,k+1))),tan(alpha)*(wakeT*b - (-1*xwHis(f,k+1))),...
                        -1*(xwHis(f,k+1)-controlPointsX(j,i)),ywHis(f,k+1)-controlPointsY(j,i),zwHis(f,k+1));

                    [u3inf,v3inf,w3inf] = VORTEX((-1*xwHis(f,k))- wakeT*b ,tan(beta)*((-1*xwHis(f,k))- wakeT*b),tan(alpha)*((-1*xwHis(f,k))- wakeT*b),...
                        wakeT*b-(-1*controlPointsX(j,i)),ywHis(f,k)+tan(beta)*(wakeT*b - (-1*xwHis(f,k)))-controlPointsY(j,i),zwHis(f,k)+tan(alpha)*(wakeT*b-(-1*xwHis(f,k))));
                    if k == 2*bPanels
                        LeadSel = m+cPanels;
                            [u6inf,v6inf,w6inf] = VORTEX(wakeT*b - (-1*xLHis(f,LeadSel)),tan(beta)*(wakeT*b - (-1*xLHis(f,LeadSel))),tan(alpha)*(wakeT*b - (-1*xLHis(f,LeadSel))),...
                            -1*(xLHis(f,LeadSel)-controlPointsX(j,i)),yLHis(f,LeadSel)-controlPointsY(j,i),zLHis(f,LeadSel));
                         u7inf =0;v7inf=0;w7inf=0;
                         u2inf=Tfac*u2inf; v2inf=Tfac*v2inf; w2inf=Tfac*w2inf;
                         u6inf=Lfac*u6inf; v6inf=Lfac*v6inf; w6inf=Lfac*w6inf;
                    elseif k == 1
                        LeadSel = m;
                             [u7inf,v7inf,w7inf] = VORTEX((-1*xLHis(f,LeadSel))- wakeT*b ,tan(beta)*((-1*xLHis(f,LeadSel))- wakeT*b),tan(alpha)*((-1*xLHis(f,LeadSel))- wakeT*b),...
                            wakeT*b-(-1*controlPointsX(j,i)),yLHis(f,LeadSel)+tan(beta)*(wakeT*b - (-1*xLHis(f,LeadSel)))-controlPointsY(j,i),zLHis(f,LeadSel)+tan(alpha)*(wakeT*b-(-1*xLHis(f,LeadSel))));
                        u6inf =0;v6inf=0;w6inf=0;
                        u3inf =Tfac*u3inf; v3inf=Tfac*v3inf; w3inf=Tfac*w3inf;
                        u7inf =Lfac*u7inf; v7inf=Lfac*v7inf; w7inf=Lfac*w7inf;
                    else
                        u6inf =0;v6inf=0;w6inf=0;u7inf =0;v7inf=0;w7inf=0;
                    end
                    end
                    end
                    %Trailing Vorticies on the wing
                    %Vortex 4: quarter chord to TE
                    %Vortex 5: TE to quarter chord
                    if k == 1
                        [u4,v4,w4] = VORTEX((-1*(panelGeomX(end)) - (-1*panelQuarterCX(m,k))),panelGeomY(end,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-controlPointsX(j,i)),panelQuarterC(m,k+1)-controlPointsY(j,i),0);
                        [u5,v5,w5] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(end))),panelQuarterC(m,k)-panelGeomY(end,k),0,-1*(panelGeomX(end)-controlPointsX(j,i)),panelGeomY(end,k)-controlPointsY(j,i),0);
                        
%                         [u4L,v4L,w4L] = VORTEX((-1*(panelGeomX(1)) - (-1*panelQuarterCX(m,k))),panelGeomY(1,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-controlPointsX(j,i)),panelQuarterC(m,k+1)-controlPointsY(i),0);
%                         [u5L,v5L,w5L] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(m))),panelQuarterC(m,k)-panelGeomY(m,k),0,-1*(panelGeomX(m)-controlPointsX(j,i)),panelGeomY(m,k)-controlPointsY(j,i),0);   
                        
%                         u5=Tfac*u5+ Lfac*u5L; v5=Tfac*v5+ Lfac*v5L; w5=Tfac*w5+ Lfac*w5L;
                         u5=Tfac*u5; v5=Tfac*v5; w5=Tfac*w5;
  
                    elseif k == 2*bPanels
%                         [u4,v4,w4] = VORTEX((-1*(panelGeomX(end)) - (-1*panelQuarterCX(m,k))),panelGeomY(end,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-controlPointsX(j,i)),panelQuarterC(m,k+1)-controlPointsY(j,i),0);
                        [u5,v5,w5] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(end))),panelQuarterC(m,k)-panelGeomY(end,k),0,-1*(panelGeomX(end)-controlPointsX(j,i)),panelGeomY(end,k)-controlPointsY(j,i),0);
                        
%                         [u4L,v4L,w4L] = VORTEX((-1*(panelGeomX(m)) - (-1*panelQuarterCX(m,k))),panelGeomY(m,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-controlPointsX(j,i)),panelQuarterC(m,k+1)-controlPointsY(j,i),0);
%                         [u5L,v5L,w5L] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(1))),panelQuarterC(m,k)-panelGeomY(1,k),0,-1*(panelGeomX(1)-controlPointsX(j,i)),panelGeomY(1,k)-controlPointsY(i),0);   
                        
%                         u4 =Tfac*u4 + Lfac*u4L; v4 =Tfac*v4 + Lfac*v4L; w4=Tfac*w4+ Lfac*w4L; 
                         u4 =0; v4 =0; w4=0; 
                    else
                        [u4,v4,w4] = VORTEX((-1*(panelGeomX(end)) - (-1*panelQuarterCX(m,k))),panelGeomY(end,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-controlPointsX(j,i)),panelQuarterC(m,k+1)-controlPointsY(j,i),0);
                        [u5,v5,w5] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(end))),panelQuarterC(m,k)-panelGeomY(end,k),0,-1*(panelGeomX(end)-controlPointsX(j,i)),panelGeomY(end,k)-controlPointsY(j,i),0);
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
gammaMatrix = reshape(gamma,cPanels,2*sum(bPanels));

% %Wake Tracking Process
% %Reset wake velocity field purturbations

%Loop through each tracked line vortex
for i = 1:2*bPanels+1
    for j=1:timeStep
        %Loop through each panel on the wing. Calculate downwash contributions BY this
        %panel 
    for k = 1:2*bPanels
        for m = 1:cPanels
                    %Bound Vortex: Vortex 1
                  [u1,v1,w1] = VORTEX(0,panelQuarterC(m,k+1)-panelQuarterC(m,k),0,-1*(panelQuarterCX(m,k)-xwHis(j,i)),panelQuarterC(m,k)-ywHis(j,i),-zwHis(j,i));
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

                   if k == 1||k==2*bPanels
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
                                -1*(xwHis(f-1,k+1)-xwHis(j,i)),ywHis(f-1,k+1)- ywHis(j,i),zwHis(f-1,k+1)-zwHis(j,i));

                            [u3,v3,w3] = VORTEX(-1*(xwHis(f-1,k)- xwHis(f,k)),ywHis(f-1,k)-ywHis(f,k),zwHis(f-1,k)-zwHis(f,k),...
                                -1*(xwHis(f,k)-xwHis(j,i)),ywHis(f,k)- ywHis(j,i),zwHis(f,k)-zwHis(j,i));

                             
                        
                            if k == 2*bPanels
                                LeadSel = m+cPanels;
                                [u6,v6,w6] = VORTEX(-1*(xLHis(f,LeadSel) - xLHis(f-1,LeadSel)),yLHis(f,LeadSel)-yLHis(f-1,LeadSel),zLHis(f,LeadSel)-zLHis(f-1,LeadSel),...
                                    -1*(xLHis(f-1,LeadSel)-xwHis(j,i)),yLHis(f-1,LeadSel)- ywHis(j,i),zLHis(f-1,LeadSel)-zwHis(j,i));
                                 u7=0; v7=0; w7=0;
                                 u6 =u6*Lfac; v6 =v6*Lfac; w6 =w6*Lfac; 
                                 u2 =u2*Tfac; v2 =v2*Tfac; w2=w2*Tfac;
                            elseif k== 1
                                LeadSel = m;
                                [u7,v7,w7] = VORTEX(-1*(xLHis(f-1,LeadSel)- xLHis(f,LeadSel)),yLHis(f-1,LeadSel)-yLHis(f,LeadSel),zLHis(f-1,LeadSel)-zLHis(f,LeadSel),...
                                    -1*(xLHis(f,LeadSel)-xwHis(j,i)),yLHis(f,LeadSel)- ywHis(j,i),zLHis(f,LeadSel)-zwHis(j,i));  
                                u6 =0; v6 =0; w6 =0; 
                                u7=u7*Lfac; v7=v7*Lfac; w7=w7*Lfac;
                                u3=u3*Tfac; v3=v3*Tfac; w3=w3*Tfac;
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
                                -1*(xwHis(f,k+1)-xwHis(j,i)),ywHis(f,k+1)-ywHis(j,i),zwHis(f,k+1)-zwHis(j,i));

                            [u3inf,v3inf,w3inf] = VORTEX((-1*xwHis(f,k))- wakeT*b ,tan(beta)*((-1*xwHis(f,k))- wakeT*b),tan(alpha)*((-1*xwHis(f,k))- wakeT*b),...
                                wakeT*b-(-1*xwHis(j,i)),ywHis(f,k)+tan(beta)*(wakeT*b-(-1*xwHis(f,k)))-ywHis(j,i),zwHis(f,k)+tan(alpha)*(wakeT*b-(-1*xwHis(f,k)))-zwHis(j,i));
                            if k == 2*bPanels
                                LeadSel = m+cPanels;
                                [u6inf,v6inf,w6inf] = VORTEX(wakeT*b - (-1*xLHis(f,LeadSel)),tan(beta)*(wakeT*b - (-1*xLHis(f,LeadSel))),tan(alpha)*(wakeT*b - (-1*xLHis(f,LeadSel))),...
                                    -1*(xLHis(f,LeadSel)-xwHis(j,i)),yLHis(f,LeadSel)-ywHis(j,i),zLHis(f,LeadSel)-zwHis(j,i));
                                
                                u7inf =0;v7inf=0;w7inf=0;
                                u2inf=Tfac*u2inf; v2inf=Tfac*v2inf; w2inf=Tfac*w2inf;
                                u6inf=Lfac*u6inf; v6inf=Lfac*v6inf; w6inf=Lfac*w6inf;
                            elseif k == 1
                                LeadSel = m;
                                [u7inf,v7inf,w7inf] = VORTEX((-1*xLHis(f,LeadSel))- wakeT*b ,tan(beta)*((-1*xLHis(f,LeadSel))- wakeT*b),tan(alpha)*((-1*xLHis(f,LeadSel))- wakeT*b),...
                                    wakeT*b-(-1*xwHis(j,i)),yLHis(f,LeadSel)+tan(beta)*(wakeT*b-(-1*xLHis(f,LeadSel)))-ywHis(j,i),zLHis(f,LeadSel)+tan(alpha)*(wakeT*b-(-1*xLHis(f,LeadSel)))-zwHis(j,i));
                                
                                u6inf =0;v6inf=0;w6inf=0;
                                u3inf =Tfac*u3inf; v3inf=Tfac*v3inf; w3inf=Tfac*w3inf;
                                u7inf =Lfac*u7inf; v7inf=Lfac*v7inf; w7inf=Lfac*w7inf; 
                            else
                                u6inf =0;v6inf=0;w6inf=0;u7inf =0;v7inf=0;w7inf=0;
                            end

                            
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
                    if k == 1
                        [u4,v4,w4] = VORTEX((-1*(panelGeomX(end)) - (-1*panelQuarterCX(m,k))),panelGeomY(end,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-xwHis(j,i)),panelQuarterC(m,k+1)-ywHis(j,i),-zwHis(j,i));
                        [u5,v5,w5] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(end))),panelQuarterC(m,k)-panelGeomY(end,k),0,-1*(panelGeomX(end)-xwHis(j,i)),panelGeomY(end,k)-ywHis(j,i),-zwHis(j,i));
                        
%                         [u4L,v4L,w4L] = VORTEX((-1*(panelGeomX(1)) - (-1*panelQuarterCX(m,k))),panelGeomY(1,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-controlPointsX(j,i)),panelQuarterC(m,k+1)-controlPointsY(i),0);
%                         [u5L,v5L,w5L] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(m))),panelQuarterC(m,k)-panelGeomY(m,k),0,-1*(panelGeomX(m)-xw(i)),panelGeomY(m,k)-yw(i),-zw(i));   
                        
%                         u5=Tfac*u5+ Lfac*u5L; v5=Tfac*v5+ Lfac*v5L; w5=Tfac*w5+ Lfac*w5L;
                         u5=Tfac*u5; v5=Tfac*v5; w5=Tfac*w5;
                        
                    elseif k == 2*bPanels
                        [u4,v4,w4] = VORTEX((-1*(panelGeomX(end)) - (-1*panelQuarterCX(m,k))),panelGeomY(end,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-xwHis(j,i)),panelQuarterC(m,k+1)-ywHis(j,i),-zwHis(j,i));
                        [u5,v5,w5] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(end))),panelQuarterC(m,k)-panelGeomY(end,k),0,-1*(panelGeomX(end)-xwHis(j,i)),panelGeomY(end,k)-ywHis(j,i),-zwHis(j,i));
                        
%                         [u4L,v4L,w4L] = VORTEX((-1*(panelGeomX(m)) - (-1*panelQuarterCX(m,k))),panelGeomY(m,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-xw(i)),panelQuarterC(m,k+1)-yw(i),-zw(i));
%                         [u5L,v5L,w5L] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(1))),panelQuarterC(m,k)-panelGeomY(1,k),0,-1*(panelGeomX(1)-controlPointsX(j,i)),panelGeomY(1,k)-controlPointsY(i),0);   
                        
%                         u4 =Tfac*u4 + Lfac*u4L; v4 =Tfac*v4 + Lfac*v4L; w4=Tfac*w4+ Lfac*w4L; 
                        u4 =Tfac*u4; v4 =Tfac*v4; w4=Tfac*w4;
                    else
                        [u4,v4,w4] = VORTEX((-1*(panelGeomX(end)) - (-1*panelQuarterCX(m,k))),panelGeomY(end,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-xwHis(j,i)),panelQuarterC(m,k+1)-ywHis(j,i),-zwHis(j,i));
                        [u5,v5,w5] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(end))),panelQuarterC(m,k)-panelGeomY(end,k),0,-1*(panelGeomX(end)-xwHis(j,i)),panelGeomY(end,k)-ywHis(j,i),-zwHis(j,i));
                    end
                    %Sum contributions from each segment in the u,v,w
                    %directions. Include factor of relaxation here

                    uw(j,i) = uw(j,i) + wRelaxS*gammaMatrix(m,k)*(u1 + u2tot + u3tot + u4 + u5 + u6tot + u7tot);
                    vw(j,i) = vw(j,i) + wRelaxS*gammaMatrix(m,k)*(v1 + v2tot + v3tot + v4 + v5 + v6tot + v7tot); 
                    ww(j,i) = ww(j,i) + wRelaxS*gammaMatrix(m,k)*(w1 + w2tot + w3tot + w4 + w5 + w6tot + w7tot);      

            end  
        end 
    end
end

%Loop through each tracked leading edge vortex
for i = 1:2*cPanels
    for j = 1:timeStep
    %Loop through each panel on the wing. Calculate downwash contributions BY this
    %panel 
        for k = 1:2*bPanels
            for m = 1:cPanels
                %Bound Vortex: Vortex 1
                  [u1,v1,w1] = VORTEX(0,panelQuarterC(m,k+1)-panelQuarterC(m,k),0,-1*(panelQuarterCX(m,k)-xLHis(j,i)),panelQuarterC(m,k)-yLHis(j,i),-zLHis(j,i));

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
                
               if k==1||k==bPanels
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
                            -1*(xwHis(f-1,k+1)-xLHis(j,i)),ywHis(f-1,k+1)- yLHis(j,i),zwHis(f-1,k+1)-zLHis(j,i));

                        [u3,v3,w3] = VORTEX(-1*(xwHis(f-1,k)- xwHis(f,k)),ywHis(f-1,k)-ywHis(f,k),zwHis(f-1,k)-zwHis(f,k),...
                            -1*(xwHis(f,k)-xLHis(j,i)),ywHis(f,k)- yLHis(j,i),zwHis(f,k)-zLHis(j,i));       

                        if k == 2*bPanels
                            LeadSel = m+cPanels;
                            [u6,v6,w6] = VORTEX(-1*(xLHis(f,LeadSel) - xLHis(f-1,LeadSel)),yLHis(f,LeadSel)-yLHis(f-1,LeadSel),zLHis(f,LeadSel)-zLHis(f-1,LeadSel),...
                                -1*(xLHis(f-1,LeadSel)-xLHis(j,i)),yLHis(f-1,LeadSel)- yLHis(j,i),zLHis(f-1,LeadSel)-zLHis(j,i));
                             u7=0; v7=0; w7=0;
                             u6 =u6*Lfac; v6 =v6*Lfac; w6 =w6*Lfac; 
                             u2 =u2*Tfac; v2 =v2*Tfac; w2=w2*Tfac;
                        elseif k == 1
                            LeadSel = m;
                            [u7,v7,w7] = VORTEX(-1*(xLHis(f-1,LeadSel)- xLHis(f,LeadSel)),yLHis(f-1,LeadSel)-yLHis(f,LeadSel),zLHis(f-1,LeadSel)-zLHis(f,LeadSel),...
                                -1*(xLHis(f,LeadSel)-xLHis(j,i)),yLHis(f,LeadSel)- yLHis(j,i),zLHis(f,LeadSel)-zLHis(j,i));  
                            u6 =0; v6 =0; w6 =0; 
                            u7=u7*Lfac; v7=v7*Lfac; w7=w7*Lfac;
                            u3=u3*Tfac; v3=v3*Tfac; w3=w3*Tfac;
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
                            -1*(xwHis(f,k+1)-xLHis(j,i)),ywHis(f,k+1)-yLHis(j,i),zwHis(f,k+1)-zLHis(j,i));

                        [u3inf,v3inf,w3inf] = VORTEX((-1*xwHis(f,k))- wakeT*b ,tan(beta)*((-1*xwHis(f,k))- wakeT*b),tan(alpha)*((-1*xwHis(f,k))- wakeT*b),...
                            wakeT*b-(-1*xLHis(j,i)),ywHis(f,k)+tan(beta)*(wakeT*b-(-1*xwHis(f,k)))-yLHis(j,i),zwHis(f,k)+tan(alpha)*(wakeT*b-(-1*xwHis(f,k)))-zLHis(j,i));
                        if k == 2*bPanels
                            LeadSel = m+cPanels;
                            [u6inf,v6inf,w6inf] = VORTEX(wakeT*b - (-1*xLHis(f,LeadSel)),tan(beta)*(wakeT*b - (-1*xLHis(f,LeadSel))),tan(alpha)*(wakeT*b - (-1*xLHis(f,LeadSel))),...
                                -1*(xLHis(f,LeadSel)-xLHis(j,i)),yLHis(f,LeadSel)-yLHis(j,i),zLHis(f,LeadSel)-zLHis(j,i));

                            u7inf =0;v7inf=0;w7inf=0;
                            u2inf=Tfac*u2inf; v2inf=Tfac*v2inf; w2inf=Tfac*w2inf;
                            u6inf=Lfac*u6inf; v6inf=Lfac*v6inf; w6inf=Lfac*w6inf;
                        elseif k == 1
                            LeadSel = m;
                            [u7inf,v7inf,w7inf] = VORTEX((-1*xLHis(f,LeadSel))- wakeT*b ,tan(beta)*((-1*xLHis(f,LeadSel))- wakeT*b),tan(alpha)*((-1*xLHis(f,LeadSel))- wakeT*b),...
                                wakeT*b-(-1*xLHis(j,i)),yLHis(f,LeadSel)+tan(beta)*(wakeT*b-(-1*xLHis(f,LeadSel)))-yLHis(j,i),zLHis(f,LeadSel)+tan(alpha)*(wakeT*b-(-1*xLHis(f,LeadSel)))-zLHis(j,i));

                            u6inf =0;v6inf=0;w6inf=0;
                            u3inf =Tfac*u3inf; v3inf=Tfac*v3inf; w3inf=Tfac*w3inf;
                            u7inf =Lfac*u7inf; v7inf=Lfac*v7inf; w7inf=Lfac*w7inf; 
                        else
                            u6inf =0;v6inf=0;w6inf=0;u7inf =0;v7inf=0;w7inf=0;
                        end
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
               if k == 1
                    [u4,v4,w4] = VORTEX((-1*(panelGeomX(end)) - (-1*panelQuarterCX(m,k))),panelGeomY(end,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-xLHis(j,i)),panelQuarterC(m,k+1)-yLHis(j,i),-zLHis(j,i));
                    [u5,v5,w5] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(end))),panelQuarterC(m,k)-panelGeomY(end,k),0,-1*(panelGeomX(end)-xLHis(j,i)),panelGeomY(end,k)-yLHis(j,i),-zLHis(j,i));

    %                         [u4L,v4L,w4L] = VORTEX((-1*(panelGeomX(1)) - (-1*panelQuarterCX(m,k))),panelGeomY(1,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-controlPointsX(j,i)),panelQuarterC(m,k+1)-controlPointsY(i),0);
%                     [u5L,v5L,w5L] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(m))),panelQuarterC(m,k)-panelGeomY(m,k),0,-1*(panelGeomX(m)-xLHis(i)),panelGeomY(m,k)-yLHis(i),-zLHis(i));   

%                     u5=Tfac*u5+ Lfac*u5L; v5=Tfac*v5+ Lfac*v5L; w5=Tfac*w5+ Lfac*w5L;
                    u5=Tfac*u5; v5=Tfac*v5; w5=Tfac*w5;
                elseif k == 2*bPanels
                    [u4,v4,w4] = VORTEX((-1*(panelGeomX(end)) - (-1*panelQuarterCX(m,k))),panelGeomY(end,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-xLHis(j,i)),panelQuarterC(m,k+1)-yLHis(j,i),-zLHis(j,i));
                    [u5,v5,w5] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(end))),panelQuarterC(m,k)-panelGeomY(end,k),0,-1*(panelGeomX(end)-xLHis(j,i)),panelGeomY(end,k)-yLHis(j,i),-zLHis(j,i));

%                     [u4L,v4L,w4L] = VORTEX((-1*(panelGeomX(m)) - (-1*panelQuarterCX(m,k))),panelGeomY(m,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-xL(i)),panelQuarterC(m,k+1)-yL(i),-zL(i));
    %                         [u5L,v5L,w5L] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(1))),panelQuarterC(m,k)-panelGeomY(1,k),0,-1*(panelGeomX(1)-controlPointsX(j,i)),panelGeomY(1,k)-controlPointsY(i),0);   

%                     u4 =Tfac*u4 + Lfac*u4L; v4 =Tfac*v4 + Lfac*v4L; w4=Tfac*w4+ Lfac*w4L; 
                    u4 =Tfac*u4; v4 =Tfac*v4; w4=Tfac*w4;
                else
                    [u4,v4,w4] = VORTEX((-1*(panelGeomX(end)) - (-1*panelQuarterCX(m,k))),panelGeomY(end,k+1)-panelQuarterC(m,k+1),0,-1*(panelQuarterCX(m,k)-xLHis(j,i)),panelQuarterC(m,k+1)-yLHis(j,i),-zLHis(j,i));
                    [u5,v5,w5] = VORTEX((-1*panelQuarterCX(m,k))- (-1*(panelGeomX(end))),panelQuarterC(m,k)-panelGeomY(end,k),0,-1*(panelGeomX(end)-xLHis(j,i)),panelGeomY(end,k)-yLHis(j,i),-zLHis(j,i));
                end
                %Sum contributions from each segment in the u,v,w
                %directions. Include factor of relaxation here
                uL(j,i) = uL(j,i) + wRelaxT*gammaMatrix(m,k)*(u1 + u2tot + u3tot + u4 + u5 + u6tot + u7tot);
                vL(j,i) = vL(j,i) + wRelaxT*gammaMatrix(m,k)*(v1 + v2tot + v3tot + v4 + v5 + v6tot + v7tot); 
                wL(j,i) = wL(j,i) + wRelaxT*gammaMatrix(m,k)*(w1 + w2tot + w3tot + w4 + w5 + w6tot + w7tot);

            end           
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
xwHis = xwHis + (xwK1);
ywHis = ywHis + (ywK1);
zwHis = zwHis + (zwK1);
xLHis = xLHis + xLK1;
yLHis = yLHis + yLK1;
zLHis = zLHis + zLK1;

xwHis = [xwHis;xwInit];
ywHis = [ywHis;ywInit];
zwHis = [zwHis;zwInit];
xLHis = [xLHis;xLInit];
yLHis = [yLHis;yLInit];
zLHis = [zLHis;zLInit];



%Solver 2: RK2
% Apply solver to track next wake location
% xw = xw + (xwK1+xwK2)/2;
% yw = yw + (ywK1+ywK2)/2;
% zw = zw + (zwK1+zwK2)/2;
%Update wake position history
% xwHis(timeStep+1,:) = xw;
% ywHis(timeStep+1,:) = yw;
% zwHis(timeStep+1,:) = zw;
% xLHis(timeStep+1,:) = xL;
% yLHis(timeStep+1,:) = yL;
% zLHis(timeStep+1,:) = zL;
% timeGlobal = timeGlobal + dt(timeStep);
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

% %Calculate the induced drag produced by the wing in potential flow with Treffz Plane Analysis
% gammaTreffz = sum(gammaMatrix);
% wiTreffz = zeros(1,length(controlPointsY));
% %Find downwash at the wake
% for i = 1:length(controlPointsY)
%                     wiTreffz(i) = sum(gammaTreffz.*(-1./(panelGeomY(1:end-1)-controlPointsY(i)) + 1./(panelGeomY(2:end)-controlPointsY(i))))/(4*pi);   
% end
% inducedDragTreffz = (rho)*sum(gammaTreffz.*wiTreffz);

%Plotting Section

%Plot Wing Descretization with horseshoe vortices and wake
wakeTPlot = wakeT*b/100;
drawSim3Conc(controlPointsY,controlPointsX,panelQuarterC,panelQuarterCX,panelGeomY,panelGeomX,alpha,beta,wakeTPlot,xwHis,ywHis,zwHis,xLHis,yLHis,zLHis,dt,bPanels,cPanels,MAC,'freeWake','NOanimate');

%Plot Panel Geometry and Control Points
%Create a Y coordinate matrix for the contour plot
bigY = zeros(cPanels,2*sum(bPanels));
for n = 1:cPanels
    bigY(n,:) = VcontrolPointsY(n,:);
end
% 
% %Plot a total lift contour plot
figure()
contourf(bigY,VcontrolPointsX,potentialLiftDist,90,'linecolor','none')
colorbar
title('Total lift distribution of planform');
xlabel('y');
ylabel('-x');
colormap jet
% 
% %Calculate spanwise lift and drag distributions by summing each column of the
% %lift and drag contours
potliftDistribution = sum(potentialLiftDist);
% 
% 
% %Divide by K to obtain the lift per unit span
for t = 1:sections
          spanwiseTotalLiftDist([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])  = ...
              potliftDistribution([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])/1;
end
% % Plot spanwise total Lift distribution
figure()
if NoDimLift == 1
    plot([-b/2 controlPointsY(end,:) b/2],[0 spanwiseTotalLiftDist/(.5*Uinf^2*rho*CAve) 0])
    title('Calculated spanwise Total Lift distribution');
    xlabel('y')
    ylabel('L''/q C_ave');
else
    plot([-b/2 controlPointsY(end,:) b/2],[0 spanwiseTotalLiftDist 0])
    title('Calculated spanwise Total Lift distribution');
    xlabel('y')
    ylabel('L''');
end
% 
% 
%Calculate and Print Coefficients
% LiftPot = sum(potliftDistribution);
% CLPot = 2*LiftPot/(rho*S*Uinf^2);
% fprintf('CL = %f\n',CLPot);
% 
% CDiTreffz = 2*inducedDragTreffz/(rho*Uinf^2*S);
% fprintf('CDi Treffz %f\n',CDiTreffz);
% 
% einv =  CLPot^2/(pi*realAR)/CDiTreffz;
% fprintf('Pot e = %f\n',einv);



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
h = 2*sqrt(spp*(spp-R1Length)*(spp-GLength)*(spp-R2Length))/GLength;

% h = norm(cross([Rx,Ry,Rz],[Gx,Gy,Gz]))/norm([Gx,Gy,Gz]);
MAC = 14.1421;
if h < (1E-3*MAC)
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