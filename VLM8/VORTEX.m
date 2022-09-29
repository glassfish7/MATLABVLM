function [u, v, w] = VORTEX(Gx, Gy, Gz, Rx, Ry, Rz)
x1n = 0;
y1n = 0;
z1n = 0;

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



R1XR2x = 1*((y-y1n)*(z-z2n)-(y-y2n)*(z-z1n));
R1XR2y = -1*((x-x1n)*(z-z2n)-(x-x2n)*(z-z1n));
R1XR2z = 1*((x-x1n)*(y-y2n)-(x-x2n)*(y-y1n));

R1XR2 = R1XR2x^2 + R1XR2y^2 + R1XR2z^2;


sig1 = ( (x2n-x1n)*(x-x1n) + (y2n-y1n)*(y-y1n) + (z2n-z1n)*(z-z1n))/...
    sqrt((x-x1n)^2 + (y-y1n)^2 + (z-z1n)^2);
sig2 = ( (x2n-x1n)*(x-x2n) + (y2n-y1n)*(y-y2n) + (z2n-z1n)*(z-z2n))/...
    sqrt((x-x2n)^2 + (y-y2n)^2 + (z-z2n)^2);


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
if R2 <= TOL2 ERROR = 1; end
if E2 <= TOL2 ERROR = 1; end
if E1 <= TOL2*R2 ERROR = 1; end

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