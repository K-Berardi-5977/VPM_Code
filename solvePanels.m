function [gamma, Vs, Cp, NumPan, Gamma, A, b, gamma_dS] = solvePanels(K, L, beta, S, U)
NumPan = length(K(:,1)); %indexing vartiable
TE_top = floor(NumPan/2);
TE_bottom = TE_top+1;
A = zeros(NumPan, NumPan); %variable to store Integral terms for linear system of equations
shift = floor(TE_top/3);
for i = 1:NumPan
    for j = 1:NumPan
        if (i == TE_top-shift && j ==TE_top || i== TE_top-shift && j==TE_bottom)
            A(i,j) = 1;

        elseif (i == j && i ~= TE_top-shift)
            A(i,j) = 0;

        elseif (i~= j && i~= TE_top-shift)
            A(i,j) = -K(i,j);

        end
        
    end
         
end

    

    

%========== Free Stream Terms ==========%
b = zeros(NumPan,1); %variable to store free stream term for linear system of equations
for n = 1:NumPan
    if n ~= TE_top-shift
        b(n) = -U*2*pi*cos(beta(n));
    else
        b(n) = 0;
    end
end

%========== Final Source Strength Calculations ============%
gamma = A\b %solving system of equations

%validation check, should be approximately zero
% sum_lambda = zeros(NumPan,1);
% for s = 1:length(gamma)
%     sum_lambda(s) = gamma(s)*S(s);
% end
% sum_lambda=sum(sum_lambda);
%========== Tangent Velocity and Pressure Coefficient

Vs = zeros(NumPan,1);
Cp = zeros(NumPan, 1); 

for i = 1:NumPan
    vortex_terms = 0;
    for j = 1:NumPan
        vortex_terms = vortex_terms + (gamma(j)/(2*pi))*(L(i,j)); %accounts for sum of all source terms in the tangent velocity
    end
    Vs(i) = U*sin(beta(i))+vortex_terms; %calculating surface velocity
    % if Vs(i)/U < 1E
    Cp(i) = 1-(Vs(i)/U)^2;

    
end
Gamma = sum(Vs);

for i = 1:length(gamma)
    gamma_dS(i) = (gamma(i)*S(i));
end

gamma_dS = sum(gamma_dS);
end

