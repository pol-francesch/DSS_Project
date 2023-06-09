function stm = stm_qns_roe_j2(t,sma0,ecc0,inc0,aop0)
    % Computes the STM for quasi-nonsingular relative orbital elements
    % including J2 perturbations from initial orbital element inputs

    % Constants
    mu = 3.986e14;          % Earth gravitational parameter [m^3/s^2]
    rE = 6378.127e3;        % Earth radius [m]
    J2 = 0.00108263;        % Earth's J2 coefficient
    
    % Auxiliary variables
    n = sqrt(mu/sma0^3);
    
    % Eccentricity-dependent parameters
    eta = sqrt(1-ecc0^2);
    kappa = 3/4*J2*rE^2*sqrt(mu)/(sma0^(7/2)*eta^4);
    E = 1+eta;
    F = 4+3*eta;
    G = 1/eta^2;
    
    % Inclination-dependent parameters
    P = 3*cos(inc0)^2-1;
    Q = 5*cos(inc0)^2-1;
    R = cos(inc0);
    S = sin(2*inc0);
    T_aux = sin(inc0)^2;
    U = sin(inc0);
    V = tan(inc0/2);
    W = cos(inc0/2)^2;

    % J2 time derivatives
    J2_constant = 3/4*J2*rE^2*sqrt(mu)/(sma0^(7/2)*eta^4);
    aopdot = J2_constant*(5*cos(inc0)^2-1);
    
    % Propagate over time
    tau = t; % time elapsed since initial
    aop_i = aop0;
    aop_f = aop_i + aopdot*tau;
    e_xi = ecc0*cos(aop_i);
    e_yi = ecc0*sin(aop_i);
    e_xf = ecc0*cos(aop_f);
    e_yf = ecc0*sin(aop_f);
    
    stm = [1, 0, 0, 0, 0, 0;
           -(3/2*n+7/2*kappa*E*P)*tau, 1, kappa*e_xi*F*G*P*tau, kappa*e_yi*F*G*P*tau, -kappa*F*S*tau, 0;
           7/2*kappa*e_yf*Q*tau, 0, cos(aopdot*tau)-4*kappa*e_xi*e_yf*G*Q*tau, -sin(aopdot*tau)-4*kappa*e_yi*e_yf*G*Q*tau, 5*kappa*e_yf*S*tau, 0;
           -7/2*kappa*e_xf*Q*tau, 0, sin(aopdot*tau)+4*kappa*e_xi*e_xf*G*Q*tau, cos(aopdot*tau)+4*kappa*e_yi*e_xf*G*Q*tau, -5*kappa*e_xf*S*tau, 0;
           0, 0, 0, 0, 1, 0;
           7/2*kappa*S*tau, 0, -4*kappa*e_xi*G*S*tau, 4*kappa*e_yi*G*S*tau, 2*kappa*T_aux*tau, 1];

end