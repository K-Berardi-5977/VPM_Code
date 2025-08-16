%Script to perform source panel method on NACA 0012 form airfoil
clc; clear;
%========== Importing NACA 0012 Airfoil Profile and Generating Panels ==========%
U = 1; %free stream velocity
c = 1; %chord length
alphad = 10;
alpha = alphad*(pi/180);


%function to load foil profile
[XB, YB, XC, YC, S, betaR, phiR, deltaD, betaD] = loadFoil2(c, alphad);


%========== Geometric Integral Terms ==========%

[K L] = Calc_Kij_Lij(XC, YC, XB, YB, phiR, S); %I=normal integral term, J=tangent integral term

%========== Linear SoE Solution, Pressure Coefficient ==========%

[gamma, V_s, Cp, NumPan, Gamma, A, b, gamma_dS] = solvePanels(K, L, betaR, S, U, betaD);

[Nx , Ny, Vxy, rp, psi, THETA, Cpxy_mesh] = vpm_plotstreamlines(XC, YC, XB, YB, phiR, S, gamma, U, alphad, Cp);





% figure; hold on; axis equal;
% 
% plot(XB,YB, 'b.', MarkerSize=7);
% plot(XC, YC, 'r*');
% plot(XB,YB,'k');
% % plot(x_c(indices), y_c(indices), 'bo', MarkerSize=7, MarkerFaceColor='c')
% title('Discretized Body Panels')
% xlabel('X')
% ylabel('Y')
% legend('Panel Bounds', 'Control Points')



XB(end) = [];
half_x = floor(NumPan/2);
figure; hold on;
axis([0 1 -2 1]);
plot(XC(1:half_x), Cp(1:half_x), 'ro');
plot(XC(half_x+1:end), Cp(half_x+1:end), 'bo'); 
plot(XC, Cp, 'k')
set(gca, 'YDir','reverse')
title(['Pressure Distribution on Airfoil Surface ($\alpha = ', num2str(alphad), ')$'], 'Interpreter','latex');
xlabel('X-Coordinate of Airfoil');
ylabel('Coefficient of Pressure (Cp)');
legend('Bottom Cp', 'Top Cp');






