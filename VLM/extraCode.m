%Code pieces not in use


%Calculate X velocity at control points

xVelMat = zeros(N,2*N);

for i = 1:length(controlPointsY)
    for j = 1:NLE
         %Loop through each panel on the wing. Calculate x velocity contributions BY this
        %panel 
        for k = 1:length(controlPointsY)
            for m = 1:N
                
                   
                   [cont1,~,~] = VortexVxyz(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),0,...
                       -1*(panelQuarterC(m,k)-controlPointsX(1,i)),panelGeomY(k)-controlPointsY(i),0); 
                   
                    [cont2,~,~] =  VortexVxyz(400*b - (-1*panelQuarterC(m,k+1)),0,0,-1*(panelQuarterC(m,k+1)-controlPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i),0);
                    
                    [cont3,~,~] =  VortexVxyz((-1*panelQuarterC(m,k))- 400*b ,0,0,400*b-(-1*controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i),0);
                   
                   xVelMat(j,i) = xVelMat(j,i) + gammaMatrix(j,i)*(cont1 + cont2 + cont3);
                
                
            end
        end
    end
end

xVelMat = xVelMat + Uinf*cos(alpha);



%2 loop method of calculating AIM 

               %Voriticty contribution if panel is on left wing
%                 if k <= N
%                     vorticityMatrix(panelCount,contCount) = ...
%                     -VXYZ(-1*(panelQuarterC(m,k)-panelQuarterC(m,k+1)),panelGeomY(k)-panelGeomY(k+1),-1*(panelQuarterC(m,k+1)-controlPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i))...
%                     - VXYZ(200*b - (-1*panelQuarterC(m,k)),0,-1*(panelQuarterC(m,k)-controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i))...
%                     - VXYZ((-1*panelQuarterC(m,k+1))- 200*b ,0,200*b-(-1*controlPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i))  ;
%                 
%                     contCount = contCount +1;
%                 end
%                 
%                  %Voriticty contribution if panel is on right wing
%                 if k > N
%                    vorticityMatrix(panelCount,contCount) = ...
%                    VXYZ(-1*(panelQuarterC(m,k+1)-panelQuarterC(m,k)),panelGeomY(k+1)-panelGeomY(k),-1*(panelQuarterC(m,k)-controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i))...    
%                     + VXYZ(200*b - (-1*panelQuarterC(m,k+1)),0,-1*(panelQuarterC(m,k+1)-controlPointsX(j,i)),panelGeomY(k+1)-controlPointsY(i))...
%                     + VXYZ((-1*panelQuarterC(m,k))- 200*b ,0,200*b-(-1*controlPointsX(j,i)),panelGeomY(k)-controlPointsY(i))  ;
%                
%                    contCount = contCount +1;
%                 end





