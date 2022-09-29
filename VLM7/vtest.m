 function [V] = vtest(Gx, Gy, Gz, Rx, Ry, Rz)
    % Compute velocity induced by unit strength vortex filament.
    % x is positive downstream, y is positive from centerline to right tip, z is upward.
    % Gx, Gy, Gz are the lengths of the vortex filament in the x, y, and z directions.
    % Rx, Ry, Rz is a vector from the vortex root to the control point: Rx = xroot-xctl
    R = [Rx,Ry,Rz]; G = [Gx,Gy,Gz]; R2 = dot(R,R); G2 = dot(G,G); eps = 1E-16;
    GR = dot(G,R); GXR = cross(G,R); GXR2 = dot(GXR,GXR); GPR = G+R; GPR2 = dot(GPR,GPR);
    Vtot = (GR/sqrt(R2)-(G2+GR)/sqrt(GPR2))/(4.0*pi*GXR2+eps); V = Vtot*GXR;
  end