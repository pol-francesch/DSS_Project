function [A_c, B_c] = get_reduced_model(oe_c)
    % Calculates the reduced model from Lippe/Steindorf
    % Inputs:
    %   roe     - state vector 
    %   oe_c    - singular orbit elements of chief [sma_c, ecc_c, inc_c, raan_c, aop_c, M_c]
    % Outputs:
    %   A_c     - plant matrix
    %   B       - control input matrix

    % Constants
    mu = 3.986e14;          % Earth gravitational parameter [m^3/s^2]
    rE = 6378.127e3;        % Earth radius [m]
    J2 = 0.00108263;        % Earth's J2 coefficient
    
    % Unwrap state
    sma_c = oe_c(1); ecc_c = oe_c(2); inc_c = oe_c(3);
    w_c   = oe_c(5); M_c  = oe_c(6);
    ex_c = ecc_c*cos(w_c);
    ey_c = ecc_c*sin(w_c);

    % Auxiliary variables
    n = sqrt(mu/sma_c^3); % mean motion
    eta = sqrt(1-ecc_c^2);
    C = sin(w_c);
    D = cos(w_c);
%     E = 1+eta;
%     F = 4+3*eta;
    G = 1/eta^2;
%     P = 3*cos(inc_c)^2-1;
    Q = 5*cos(inc_c)^2-1;
    %R = cos(inc_c);
    S = sin(2*inc_c);
    T = sin(inc_c)^2;
    %U = sin(inc_c);
    kappa = 3/4*J2*rE^2*sqrt(mu)/(sma_c^(7/2)*eta^4);

    % Linear dynamic model
    % Plant matrix, A_c = A_kep + A_J2 [Koenig]
    % Accounts for Keplerian motion and J2 perturbations
%     A_kep = [ [0; -1.5*n; zeros(4,1)], ...
%               [zeros(2,5);zeros(4,5)] ]; % [Koenig et al, New STMs]
    
%     A_J2 = kappa*[0, 0, 0, 0, 0, 0;
%             -7/2*eta*P, 0, 3*ecc_c*G*P, 0, -3*eta*S, 0;
%             0, 0, 0, 0, 0, 0;
%             -7/2*Q, 0, 4*ecc_c*G*Q, 0, -5*S, 0;
%             0, 0, 0, 0, 0, 0;
%             7*R, 0, -8*ecc_c*G*R, 0, 2*U, 0];

%     A_J2 = [0, 0, 0, 0, 0, 0;
%             -7/2*E*P, 0, ecc_c*F*G*P, 0, -F*S, 0;
%             0, 0, 0, 0, 0, 0;
%             -7/2*ecc_c*Q, 0, 4*ecc_c^2*G*Q, 0, -5*ecc_c *S, 0;
%             0, 0, 0, 0, 0, 0;
%             7/2*S, 0, -4*ecc_c*G*S, 0, 2*T, 0];
%     A_J2 = [[A_J2(1,1),A_J2(1,3:6)]; [A_J2((3:6),1),A_J2((3:6),3:6)] ];
%     
%     A_c = A_J2;

    % [Steindorf]
    A_c = kappa*[0, 0, 0, 0, 0;
                 7/2*ey_c*Q, -(4*ex_c*ey_c*G+C)*Q, -(1+4*ey_c^2*G-D)*Q, 5*ey_c*S, 0;
                 -7/2*ex_c*Q, (1+4*ex_c^2*G-D)*Q, (4*ex_c*ey_c*G-C)*Q, -5*ex_c*S, 0;
                 0, 0, 0, 0, 0;
                 7/2*S, -4*ex_c*G*S, -4*ey_c*G*S, 2*T, 0];

    % Control input matrix -- only use 2nd and 3rd columns
    % Near-circular [Lippe pg. 56]
%     B_circ = 1/(sma_c*n) * [0, 2, 0;
%                        -2, 0, 0;
%                        sin(u_M), 2*cos(u_M), 0;
%                        -cos(u_M), 2*sin(u_M), 0;
%                        0, 0, cos(u_M);
%                        0, 0, sin(u_M)];
    % B = B_circ;
    % B = B(:,2:3);
    % B = [B(1,:); B(3:end,:)];

    % For eccentric orbits [Steindorf pg. 5]
    u_t = mean2true(M_c, ecc_c, 1e-12); % true anomaly = mean anomaly for small eccentricity
  
    B_ecc = 1/(sma_c*n)*[2/eta*(1+ecc_c*cos(u_t)), 0;
                           eta*( (2+ecc_c*cos(u_t))*cos(w_c+u_t)+ex_c ) / (1+ecc_c*cos(u_t)),  (eta*ey_c)/tan(inc_c)*(sin(w_c+u_t)) / (1+ecc_c*cos(u_t));
                           eta*( (2+ecc_c*cos(u_t))*sin(w_c+u_t)+ey_c ) / (1+ecc_c*cos(u_t)), (-eta*ex_c)/tan(inc_c)*(sin(w_c+u_t)) / (1+ecc_c*cos(u_t));
                           0, eta*(cos(w_c+u_t)) / (1+ecc_c*cos(u_t));
                           0, eta*(sin(w_c+u_t)) / (1+ecc_c*cos(u_t))];
    B_c = B_ecc;

end