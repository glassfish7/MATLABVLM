function [u,v,w] = VORTEXVec2(Gx, Gy, Gz, Rx, Ry, Rz)
global timeGlobal
global MAC
global spanGlobal
%Vortex Root at Origin
x1n = 0;
y1n = 0;
z1n = 0;

%Vortex Tip Location
x2n = Gx;
y2n = Gy;
z2n = Gz;

% Set Control Point location

% Convention 1: Textbook convention(Control Point minus Vortex Root)
% x = Rx;
% y = Ry;
% z = Rz;

%Convention 2: Kroo convention(Vortex Root minus Control Point)
x = -Rx;
y = -Ry;
z = -Rz;
%R1 X R2
R1XR2x = 1*((y-y1n).*(z-z2n)-(y-y2n).*(z-z1n));
R1XR2y = -1*((x-x1n).*(z-z2n)-(x-x2n).*(z-z1n));
R1XR2z = 1*((x-x1n).*(y-y2n)-(x-x2n).*(y-y1n));

%|R1 X R2|
R1XR2 = R1XR2x.^2 + R1XR2y.^2 + R1XR2z.^2;

%sigma: see Virginia Tech textbook ch 6 pg. 18
sig1 = ( (x2n-x1n).*(x-x1n) + (y2n-y1n).*(y-y1n) + (z2n-z1n).*(z-z1n))./...
    sqrt((x-x1n).^2 + (y-y1n).^2 + (z-z1n).^2);
sig2 = ( (x2n-x1n).*(x-x2n) + (y2n-y1n).*(y-y2n) + (z2n-z1n).*(z-z2n))./...
    sqrt((x-x2n).^2 + (y-y2n).^2 + (z-z2n).^2);

%Apply Finite Core Model
G2=Gx.*Gx+Gy.*Gy+Gz.*Gz;
R2=Rx.*Rx+Ry.*Ry+Rz.*Rz;
E2=(Gx+Rx).*(Gx+Rx)+(Gy+Ry).*(Gy+Ry)+(Gz+Rz).*(Gz+Rz);
R1Length = sqrt(R2);
GLength = sqrt(G2);
R2Length = sqrt(E2);
spp = (R1Length+GLength+R2Length)/2;
h = 2*GLength.*sqrt(spp.*(spp-R1Length).*(spp-GLength).*(spp-R2Length));
% h = norm(cross([Rx,Ry,Rz],[Gx,Gy,Gz]))/norm([Gx,Gy,Gz]);


V = (sig1-sig2) / (4.0*pi);
u = V.*(R1XR2x./R1XR2);
v = V.*(R1XR2y./R1XR2);
w = V.*(R1XR2z./R1XR2);


% h < sqrt(4*.001*timeGlobal)
u(h < 1E-3*spanGlobal) = 0;
v(h < 1E-3*spanGlobal) = 0;
w(h < 1E-3*spanGlobal) = 0;




%Apply ArtViscCore
% G2=Gx*Gx+Gy*Gy+Gz*Gz;
% R2=Rx*Rx+Ry*Ry+Rz*Rz;
% E2=(Gx+Rx)*(Gx+Rx)+(Gy+Ry)*(Gy+Ry)+(Gz+Rz)*(Gz+Rz);
% R1Length = sqrt(R2);
% GLength = sqrt(G2);
% R2Length = sqrt(E2);
% spp = (R1Length+GLength+R2Length)/2;
% h = 2*GLength*sqrt(spp*(spp-R1Length)*(spp-GLength)*(spp-R2Length));
% MAC = 9.4412;
% 
% V = (sig1-sig2) / (4.0*pi);
% % disp(timeGlobal)
% V = V*(1-exp(-(h)^2/(4*(.001*MAC)*(timeGlobal))));
% u = V*(R1XR2x/R1XR2);
% v = V*(R1XR2y/R1XR2);
% w = V*(R1XR2z/R1XR2);


end