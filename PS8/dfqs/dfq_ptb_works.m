function zdot = dfq_ptb_works(t,z,r_main,mu,J2,u_TN, eci_d)
    % Calculates the derivative of z at t with J2 (FODE)
    % Inputs:
    %   t       - current time, s
    %   z       - state vector [rx, ry, rz, vx, vy, vz] [m, m/s]
    %   r_main  - main body radius, m
    %   mu      - gravitaional parameter of central body, m^3/s^2
    %   J2      - J2 coefficient
    %   u_TN    - TN control inputs
    %   oe_c    - deputy orbit elements [sma_c, ecc_c, inc_c, raan_c, aop_c, M_c]
    % Outputs:
    %   zdot    - derivative of state vector [vx, vy, vz, ax, ay, az]

    % current state quantities
    r = z(1:3);
    rdot = z(4:6);
    r_mag = norm(r);

    % J2 propagation
    x = r(1); y = r(2); z = r(3);    
    J2_perturb = 3*J2*mu*r_main^2/(2*r_mag^5) * [(5*z^2/r_mag^2-1)*x;
                                                 (5*z^2/r_mag^2-1)*y;
                                                 (5*z^2/r_mag^2-3)*z];
    
    % Control
    u_RTN = [0,u_TN]';     % Full RTN
    eci2rtn = get_hill_frame(eci_d);
    u_ECI = (eci2rtn') * u_RTN;

    % Total perturbation
    d = J2_perturb + u_ECI;
    
    % EOM
    rddot = -mu*r/r_mag^3 + d;
    zdot = [rdot; rddot];
end