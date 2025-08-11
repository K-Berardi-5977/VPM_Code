function [Px, Py, dx, dy, S, phi, deltaD, beta, betaD, phiD] = loadFoil(c, t_max, alphad)

load("c_x.mat"); %profile x/c coordinates
load("t_x.mat"); %profile
Px = c.*c_x; %for bottom surface 
Pxb = flip(Px); %creates second half of curve s (x-coordinate)
Px = vertcat(Px, Pxb); %for bottom surface
Py = t_max*c.*t_x; %ypoints generated based on max thickness ratio to chord length
Py = Py/2; %halfing y-values since they are based on thickness (i.e., "diameter") 
Pyb = -1*flip(Py);
Py = vertcat(Py, Pyb);


dx = diff(Px);
dy = diff(Py);


for i = 1:length(dx)
    S(i) = sqrt(dx(i)^2 + dy(i)^2);
    phiD(i)=atan2d(dy(i), dx(i));
    if phiD(i) < 0
        phiD(i) = phiD(i) + 360;
    end


end
deltaD = phiD + 90;
betaD = deltaD - alphad;
phi = phiD.*(pi/180);
betaD(betaD>360) = betaD(betaD>360) - 360;
beta = betaD.*(pi/180);
end