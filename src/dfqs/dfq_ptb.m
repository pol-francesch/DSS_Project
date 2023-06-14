function zdot = dfq_ptb(t,z,r_main,mu,J2,u_TN, eci_d)
    % Calculates the derivative of z at t with J2 (FODE)
    % Inputs:
    %   t       - current time, s
    %   z       - state vector [rx, ry, rz, vx, vy, vz] [m, m/s]
    %   r_main  - main body radius, m
    %   mu      - gravitaional parameter of central body, m^3/s^2
    %   J2      - J2 coefficient
    %   u       - RTN control inputs (accelerations)
    %   eci_d   - deputy position/velocity in ECI [rx, ry, rz, vx, vy, vz]
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
    u_RTN = u_TN';     % Full RTN
    eci2rtn = get_hill_frame(eci_d);
    u_ECI = (eci2rtn') * u_RTN;

    % Total perturbation
    d = J2_perturb + u_ECI;
    
    % EOM
    rddot = -mu*r/r_mag^3 + d;
    zdot = [rdot; rddot];
end