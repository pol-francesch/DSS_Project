function [dV, error] = least_squares_maneuvers2(mu, J2, rE, N, u_ks, delta_us, desired_delta_roe, abs_chief_oe)
    % Given maneuver locations, solve least squares to get the maneuver
    % delta v and the error

    sma_c = abs_chief_oe(1);
    ecc_c = abs_chief_oe(2);
    inc_c = abs_chief_oe(3);
    aop_c = abs_chief_oe(5);

    n = sqrt(mu/sma_c^3);
    T = 2*pi*sqrt(sma_c^3/mu); 

    % Control input matrices
    gamma = @(u_k) 1/(n*sma_c)*[0, 2, 0;
                             -2, 0, 0;
                             sin(u_k), 2*cos(u_k), 0;
                             -cos(u_k), 2*sin(u_k), 0;
                             0, 0, cos(u_k);
                             0, 0, sin(u_k)];
%     % Simplified dynamics
    eta = sqrt(1-ecc_c^2);
    kappa = 3/4*J2*rE^2*sqrt(mu)/(sma_c^(7/2)*eta^4);
    P = 3*cos(inc_c)^2-1;
    Q = 5*cos(inc_c)^2-1;
%     E = 1+eta;
%     F = 4+3*eta;
%     S = sin(2*inc_c);
%     J2_constant = 3/4*J2*rE^2*sqrt(mu)/(sma_c^(7/2)*eta^4);
%     aopdot = J2_constant*(5*cos(inc_c)^2-1);
    tau = @(delta_u) delta_u/(n+kappa*(eta*P+Q));
%     
%     % Specialized stm for near-circular chief orbit
%     stm = @ (tau_k) [1, 0, 0, 0, 0, 0;
%             -(7*kappa*E*P+3*n)/2*tau_k, 1, 0, 0, -kappa*F*S*tau_k, 0;
%             0, 0, cos(aopdot*tau_k), -sin(aopdot*tau_k), 0, 0;
%             0, 0, sin(aopdot*tau_k), cos(aopdot*tau_k), 0, 0;
%             0, 0, 0, 0, 1, 0;
%             7/2*kappa*S*tau_k, 0, 0, 0, 2*kappa*T*tau_k, 1];

    % Build maneuvers matrix
    M = [];
    
    for i=1:N
        M = [M, stm2(tau(delta_us(i)),J2,rE,mu,sma_c,ecc_c,inc_c,aop_c) * gamma(u_ks(i))];
    end
    
    % Least squares to get the dv
    [dV,~] = lsqr(M, desired_delta_roe);
    
    droe_actual = M*dV;
    error = droe_actual - desired_delta_roe;
end


function stm = stm2(tau,J2,rE,mu,sma_c,ecc_c,inc_c,aop_c)

    % Eccentricity-dependent parameters
    eta = sqrt(1-ecc_c^2);
    n = sqrt(mu/sma_c^3);
    kappa = 3/4*J2*rE^2*sqrt(mu)/(sma_c^(7/2)*eta^4);
    E = 1+eta;
    F = 4+3*eta;
    G = 1/eta^2;
    
    % Inclination-dependent parameters
    P = 3*cos(inc_c)^2-1;
    Q = 5*cos(inc_c)^2-1;
    R = cos(inc_c);
    S = sin(2*inc_c);
    T_aux = sin(inc_c)^2;
    U = sin(inc_c);
    V = tan(inc_c/2);
    W = cos(inc_c/2)^2;

    % J2 time derivatives
    J2_constant = 3/4*J2*rE^2*sqrt(mu)/(sma_c^(7/2)*eta^4);
    aopdot = J2_constant*(5*cos(inc_c)^2-1);
    aop_i = aop_c;
    aop_f = aop_i + aopdot*tau;
    e_xi = ecc_c*cos(aop_i);
    e_yi = ecc_c*sin(aop_i);
    e_xf = ecc_c*cos(aop_f);
    e_yf = ecc_c*sin(aop_f);

    stm = [1, 0, 0, 0, 0, 0;
           -(3/2*n+7/2*kappa*E*P)*tau, 1, kappa*e_xi*F*G*P*tau, kappa*e_yi*F*G*P*tau, -kappa*F*S*tau, 0;
           7/2*kappa*e_yf*Q*tau, 0, cos(aopdot*tau)-4*kappa*e_xi*e_yf*G*Q*tau, -sin(aopdot*tau)-4*kappa*e_yi*e_yf*G*Q*tau, 5*kappa*e_yf*S*tau, 0;
           -7/2*kappa*e_xf*Q*tau, 0, sin(aopdot*tau)+4*kappa*e_xi*e_xf*G*Q*tau, cos(aopdot*tau)+4*kappa*e_yi*e_xf*G*Q*tau, -5*kappa*e_xf*S*tau, 0;
           0, 0, 0, 0, 1, 0;
           7/2*kappa*S*tau, 0, -4*kappa*e_xi*G*S*tau, 4*kappa*e_yi*G*S*tau, 2*kappa*T_aux*tau, 1];
end
