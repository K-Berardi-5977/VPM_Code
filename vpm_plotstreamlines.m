function [Nx , Ny, Vxy, rp, psi, THETA] = vpm_plotstreamlines(xi, yi, Xj, Yj, phi, S, gamma, U, alphad, Cp)
numPan = length(xi);
[Nx, Ny] = calc_NxNy(xi, yi, Xj, Yj, phi, S)

nGridX = numPan;                                                           % X-grid for streamlines and contours
nGridY = numPan;                                                           % Y-grid for streamlines and contours
xVals  = [-0.3; 1.3];                                                   % X-grid extents [min, max]
yVals  = [-.8; 0.8];                                                   % Y-grid extents [min, max]

% Generate the grid points
Xgrid   = linspace(xVals(1),xVals(2),nGridX)';                          % X-values in evenly spaced grid
Ygrid   = linspace(yVals(1),yVals(2),nGridY)';                          % Y-values in evenly spaced grid
[XX,YY] = meshgrid(Xgrid,Ygrid);                                        % Create meshgrid from X and Y grid arrays
rr = sqrt(XX.^2 + YY.^2);
THETA = atan2d(YY, XX);
% Streamline parameters
stepsize = 0.01;                                                        % Step size for streamline propagation
maxVert  = nGridX*nGridY*10;                                            % Maximum vertices
slPct    = 30;                                                          % Percentage of streamlines of the grid
Ysl      = linspace(yVals(1),yVals(2),floor((slPct/100)*nGridY))';      % Create array of Y streamline starting points



% Initialize velocities
Vx = zeros(nGridX,nGridY);                                              % Initialize X velocity matrix
Vy = zeros(nGridX,nGridY);

for m = 1:1:nGridX   %iterating over the jth control point
    for n = 1:1:nGridY  %for each control point, iterate over j=1:n panels
            XP = XX(m,n);   %Current iteration's X grid point
            YP = YY(m,n);   %Current iteration's Y grid point
            [Nxx, Nyy] = calc_NxNy(XP, YP, Xj, Yj, phi,S)
            [in,on] = inpolygon(XP,YP,Xj, Yj);
             % See if points are in or on the airfoil
             if (in == 1 || on == 1)                                         % If the grid point is in or on the airfoil
                 Vx(m,n) = 0;                                                % Set X-velocity equal to zero
                 Vy(m,n) = 0;                                                % Set Y-velocity equal to zero
            else                                                            % If the grid point is outside the airfoil
                Vx(m,n) = U*cosd(alphad) + sum((gamma.*Nxx)./(2*pi));         % Compute X-velocity
                Vy(m,n) = U*sind(alphad) + sum((gamma.*Nyy)./(2*pi));         % Compute Y-velocity
            end
          
    end

end
Vxy =  sqrt(Vx.^2 + Vy.^2);
gamma_dS = gamma(:).*S(:);
psi = U*(YY*cosd(alphad)-XX*sind(alphad));

for j=1:length(gamma_dS)
    dxp = (XX-xi(j));
    dyp = (YY-yi(j));
    rp = sqrt(dxp.^2 + dyp.^2);
    psi = psi + (gamma_dS(j)/(2*pi))*log(rp);
end

[row, col] = find(Cp>0.99945)

in = inpolygon(XX, YY, Xj, Yj);
psi_mask = psi;
psi_mask(in) = NaN;

figure; hold on; box on
contour(XX, YY, psi_mask, 60, 'LineWidth', 0.9);
fill(Xj, Yj, [0.15 0.15 0.15], 'EdgeColor', 'k'); % the body
for i = 1:length(row)
    plot(Xj(row(i)), Yj(row(i)), 'ro', 'MarkerSize', 3, MarkerFaceColor='r')
end
axis equal tight
xlabel('x/c'); ylabel('y/c');
title(['VPM Streamlines ($\alpha = ', num2str(alphad), '$ deg)' ], 'Interpreter','latex')
legend("Streamlines",'Airfoil Body','Stagnation Points')
end

