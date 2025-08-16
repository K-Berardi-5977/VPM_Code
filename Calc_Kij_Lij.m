function [K, L] = Calc_Kij_Lij(xi, yi, Xj, Yj, phi, S)
NumPan = length(xi); %iteration variable
K = zeros(NumPan, NumPan);
L = zeros(NumPan, NumPan);

% xi & yi - vectors containing the x and y coordinates of the control
% points

%Defining Integral Terms:
for i = 1:NumPan %iterating over the ith control point
    for j = 1:NumPan %for each control point, iterate over j=1:n panels
        if j ~= i
            %all of the coefficient terms determined in the solution of the
            %geometric integrals 
            A = -(xi(i)-Xj(j))*cos(phi(j))-(yi(i)-Yj(j))*sin(phi(j));
            B = (xi(i)-Xj(j))^2 +(yi(i)-Yj(j))^2;
            Cn = -cos(phi(i)-phi(j)); %C coefficient for normal velocity integral term
            Dn = (xi(i)-Xj(j))*cos(phi(i)) + (yi(i)-Yj(j))*sin(phi(i)); %D coefficient for normal velocity integral term
            Cs = sin(phi(j)-phi(i)); %C coefficient for tangent(surface) velocity integral term
            Ds = (xi(i)-Xj(j))*sin(phi(i))-(yi(i)-Yj(j))*cos(phi(i)); %D coefficient for tangent(surface) velocity integral term
            E = sqrt(B-A^2);
            Sj = S(j);

            if ~isreal(E)
                E = 0;
            end
            %compute normal geometric integral for each panel at the ith
            %point
            K(i,j)= (Cn/2)*log((Sj^2 + 2*A*Sj + B)/B) + ((Dn - A*Cn)/E)*(atan2(Sj+A, E) ...
                -atan2(A, E)); 

            %compute tangent (surface) geometric interal for each panel at
            %the ith point
            L(i,j) = (Cs/2)*log((Sj^2 + 2*A*Sj + B)/B) + ((Ds - A*Cs)/E)*(atan2(Sj+A, E) ...
                -atan2(A, E));

       
    if (isnan(K(i,j)) || isinf(K(i,j)) || ~isreal(K(i,j)))
        K(i,j) = 0;
    end

 

    if (isnan(L(i,j)) || isinf(L(i,j)) || ~isreal(L(i,j)))
        L(i,j) = 0;
    end

    end
   
end
end
end