function [gamma, Vs, Cp, NumPan, Gamma, A, b, gamma_dS] = solvePanels(K, L, beta, S, U)

NumPan = length(K(:,1)); %indexing vartiable


%========== Vortex Strength Influence Coefficients ==========%
A = zeros(NumPan, NumPan); %variable to store influence coefficients for linear system of equations

for i = 1:NumPan
    for j = 1:NumPan
        if (i == 20&& j == 1) || (i ==20 && j == NumPan) %enforce the Kutta condition at selected panel 
            A(i,j) = 1;
        elseif (i == j && i ~= 20)
            A(i,j) = 0; %no self influence 
        elseif (i~= j && i~= 20)
            A(i,j) = -K(i,j);
        end
    end 
end


%========== Free Stream Terms ==========%
b = zeros(NumPan,1); %variable to store free stream term for linear system of equations
for n = 1:NumPan
    if n ~= 20
        b(n) = -U*2*pi*cos(beta(n));
    else
        b(n) = 0;
    end
end

%========== Final Source Strength Calculations ============%
gamma = A\b %solving system of equations

%========== Surface Velocity and Pressure Coefficient ==========%
Vs = zeros(NumPan,1); %initialize surface velocity vector
Cp = zeros(NumPan, 1); %initialize pressure coefficient vector

for i = 1:NumPan
    vortex_terms = 0;
    for j = 1:NumPan
        vortex_terms = vortex_terms + (gamma(j)/(2*pi))*(L(i,j)); %accounts for sum of all source terms in the tangent velocity
    end
    Vs(i) = U*sin(beta(i))-vortex_terms; %calculating surface velocity at ith panel control point
    Cp(i) = 1-(Vs(i)/U)^2; %calculating pressure coefficient at ith panel control point
end

Gamma = sum(Vs); %validation check - from definition of circulation

%calculate the strength of each panel
for i = 1:NumPan
    gamma_dS(i) = (gamma(i)*S(i));
end

gamma_dS = sum(gamma_dS);  %validation check = 'Gamma' where 'Gamma' is the sum of the flow velocity along the surface (circulation)
end

