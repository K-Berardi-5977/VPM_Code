function [XB, YB, XC, YC, S, betaR, phiR, deltaD] = loadFoil2(c, alphad)


foilShape = load('naca0012.dat');

XB = foilShape(:,1); %profile X-coordinates
YB = foilShape(:,2); %profile Y-coordinates

XB = flip(XB);
YB = flip(YB);
% np = 67; %number of boundary points
% 
% XB_top = linspace(0,1,np)';
% XB_bottom = flip(XB_top);
% XB_bottom(1) = [];
% XB = vertcat(XB_bottom, XB_top)
% 
% for i = 1:np
%     YB_top(i,1) = 0.594689181*(0.298222773*sqrt(XB(i)) - 0.127125232*XB(i) - 0.357907906*XB(i)^2 + 0.291984971*XB(i)^3 - 0.105174606*XB(i)^4)
% 
% end
% YB_bottom = -YB_top
% YB_bottom = flip(YB_bottom);
% YB_bottom(1) = [];
% YB = vertcat(YB_top, YB_bottom)
numPan = length(YB)-1; %N-1 panels for N coordinate points



edge = zeros(numPan, 1);
for i = 1:numPan
    edge(i) = (XB(i+1)-XB(i))*(YB(i+1)+YB(i)); %check for direction of boundary points
end
sumEdge = sum(edge);

% flip array if panels are oriented CCW
if (sumEdge < 0)
    flipud(XB);
    flipud(YB);
end

%initialize control point, panel length, and orientation variables 
XC = zeros(numPan, 1); %X control point coordinate
YC = zeros(numPan, 1); %Y control point coordinate
S = zeros(numPan, 1); %panel lengths
phiD = zeros(numPan, 1); %panel orientation

dX = zeros(numPan, 1); %XB(i+1)-XB(i) --- for use in panel orientation and length
dY = zeros(numPan, 1); %YB(i+1)-YB(i) --- for use in panel orientation and length

for i = 1:numPan
    XC(i) = 0.5*(XB(i+1)+XB(i)); %X-coordinate of control point of ith panel
    YC(i) = 0.5*(YB(i+1)+YB(i)); %Y-coordinate of control point of ith panel
    dX(i) = (XB(i+1)-XB(i)); %change in x along ith panel
    dY(i) = (YB(i+1)-YB(i)); %change in y along ith panel
    S(i) = sqrt(dX(i)^2 + dY(i)^2); %length of ith panel
    phiD(i) = atan2d(dY(i), dX(i)); %angle of ith panel with respect to positive x-axis
    if (phiD(i) < 360)
        phiD(i) = phiD(i) + 360; %makes all angles positive angles 
    end
end

deltaD = phiD + 90; %angle from positive x-axis to outward facing unit normal 
betaD = deltaD-alphad; %angle between freestream velocity vector and outward facing unit normal 
betaD(betaD>360) = betaD(betaD>360)-360 %make angles less than 360

phiR = phiD.*(pi/180); %angle phi converted to radians
betaR = betaD.*(pi/180); %angle beta converted to radians


end