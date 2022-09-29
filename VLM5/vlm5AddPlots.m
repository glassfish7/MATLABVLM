%VLM 5.0: Additional Plotting Options
%Run this after running vlm5 to get additional plots
%All data from vlm5 must be in memory
%Plot a potential lift contour plot
figure()
contourf(bigY,VcontrolPointsX,potentialLiftDist,90,'linecolor','none')
colorbar
title('Potential lift contour plot of planform');
xlabel('y');
ylabel('-x');
colormap jet


%Plot a vortex lift contour plot
figure()
contourf(bigY,VcontrolPointsX,vortexLift,90,'linecolor','none')
colorbar
title('Vortex lift distribution of planform');
xlabel('y');
ylabel('-x');
colormap jet



%Calculate spanwise potential lift distribution by summing each column of the
%lift distribution
potliftDistribution = sum(potentialLiftDist);
NVLiftDistribution = sum(combineLiftNV);

%Divide by K to obtain the lift per unit span
for t = 1:sections
          spanwiseLiftDist([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])  = ...
              potliftDistribution([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])/K(t);
       
          spanwiseNVLiftDistribution([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])  = ...
              NVLiftDistribution([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])/K(t);
          
          spanwisePotIndDragDistribution([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])  = ...
              potIndDragDistribution([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])/K(t);
          
          spanwiseNVDist([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])  = ...
              NVDist([sPoints(1,t):ePoints(1,t),sPoints(2,t):ePoints(2,t)])/K(t);
          
end

% Plot spanwise lift distribution
figure()
if NoDimLift == 1
    plot([-b/2 controlPointsY b/2],[0 spanwiseLiftDist/(.5*Uinf^2*rho*CAve) 0])
    title('Calculated spanwise Potential Flow Lift distribution');
    xlabel('y')
    ylabel('L''/q C_ave');
else
    plot([-b/2 controlPointsY b/2],[0 spanwiseLiftDist 0])
    title('Calculated spanwise Potential Flow Lift distribution');
    xlabel('y')
    ylabel('L''');
end



% Plot spanwise no suction lift distribution
figure()
if NoDimLift == 1
    plot([-b/2 controlPointsY b/2],[0 spanwiseNVLiftDistribution/(.5*Uinf^2*rho*CAve) 0])
    title('Calculated spanwise zero suction zero vortex Lift Distribution');
    xlabel('y')
    ylabel('L''/q C_ave');
else
    plot([-b/2 controlPointsY b/2],[0 spanwiseNVLiftDistribution 0])
    title('Calculated spanwise zero suction zero vortex Lift Distribution');
    xlabel('y')
    ylabel('L''');
end

%Plot spanwise potential drag
figure()
plot([-b/2 controlPointsY b/2],[0 spanwisePotIndDragDistribution 0])
title('Calculated spanwise potential induced drag distribution');
xlabel('y')
ylabel('D''');

%Plot spanwise no suction drag
figure()
plot([-b/2 controlPointsY b/2],[0 spanwiseNVDist 0])
title('Calculated spanwise no suction induced drag distribution');
xlabel('y')
ylabel('D''');
