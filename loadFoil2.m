function [XB, YB, XC, YC, S, betaR, phiR, deltaD] = loadFoil2(c, t_max, alphad)

load("c_x_VPM.mat"); %profile x/c coordinates
load("t_x_VPM.mat"); %profile

numB = length(t_x); %iteration variable for number of boundarypoints

XB = zeros(numB,1); %initialize x-coordinate of boundary points
YB = zeros(numB,1); %initialize y-coordinate of boundary points

for i = 1:numB
    XB(i) = c*c_x(i);
    YB(i) = t_max*c*t_x(i);
end

X_te_panels1 = linspace(XB(end-1) ,XB(end), 20)';
Y_te_panels1 = linspace(YB(end-1), YB(end), 20)';

% X_te_panels1(1) = [];
% Y_te_panels1(1) = [];
% 
% X_te_panels2 = linspace(XB(end-2), XB(end-1), 10)';
% Y_te_panels2 = linspace(YB(end-2),YB(end-1), 10)';


XB = vertcat(  XB(1:end-2), X_te_panels1)
YB = vertcat( YB(1:end-2), Y_te_panels1)

XB_bottom = XB(1:end-1);
size(XB)
size(XB_bottom)
XB_bottom = flip(XB_bottom);
XB = vertcat(XB, XB_bottom);


YB_bottom = YB(1:end-1);
YB_bottom = -1.*(flip(YB_bottom));
YB = vertcat(YB, YB_bottom)


numPan = length(YB) - 1; %number of panels



edge = zeros(numPan, 1);
for i = 1:numPan
    edge(i) = (XB(i+1)-XB(i))*(YB(i+1)+YB(i)); %checl for direction of boundary points
end
sumEdge = sum(edge);
%
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