function [Vx, Vy, Vz] = VortexVxyz(Gx, Gy, Gz, Rx, Ry, Rz)
    % Compute velocity induced by a unit strength vortex filament.
    % Vx, Vy, Vz are the velocity components in the coordinate directions.
    % Gx, Gy, Gz are the vortex lengths in each direction: tip - root.
    % Rx, Ry, Rz is the radius vector from vortex root to control point.
    

    TOL = 1.0E-10;
    R2=Rx*Rx+Ry*Ry+Rz*Rz;
    G2=Gx*Gx+Gy*Gy+Gz*Gz;
    GR=Gx*Rx+Gy*Ry+Gz*Rz;
    TOL2 = TOL * G2;
 

    GXRx=Gy*Rz-Gz*Ry;
    GXRy=Gz*Rx-Gx*Rz;
    GXRz=Gx*Ry-Gy*Rx;
 

    % Check to see if control point lies in line with vortex:
    E1=GXRx*GXRx+GXRy*GXRy+GXRz*GXRz;
    E2=(Gx+Rx)*(Gx+Rx)+(Gy+Ry)*(Gy+Ry)+(Gz+Rz)*(Gz+Rz);
    ERROR = 0;
    if R2 <= TOL2 ERROR = 1; end
    if E2 <= TOL2 ERROR = 1; end
    if E1 <= TOL2*R2 ERROR = 1; end
 

    if ERROR==0
       V = (GR/sqrt(R2) - (G2+GR)/sqrt(E2)) / (4.0*pi*E1);
       Vx = V*GXRx;
       Vy = V*GXRy;
       Vz = V*GXRz;
    else
       Vx=0.0;
       Vy=0.0;
       Vz=0.0;
    end
 

end