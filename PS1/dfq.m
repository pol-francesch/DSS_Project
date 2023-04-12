function zdot = dfq(t,z, r_main, mu, J2)
    % Calculates the derivative of z at t wih J2
    % Inputs:
    %   t       - current time, s
    %   z       - state vector [rx, ry, rz, vx, vy, vz] [m, m/s]
    %   r_main  - main body radius, m
    %   mu      - gravitaional parameter of central body, m^3/s^2
    %   J2      - J2 coefficient
    % Outputs:
    %   zdot    - derivative of state vector [vx, vy, vz, ax, ay, az]

    % current state quantities
    r = z(1:3);
    rdot = z(4:6);
    r_mag = norm(r);

    % J2 propagation
    x = r(1); y = r(2); z = r(3);
    % J2_perturb = [-mu*x/r_mag^3*J2*3/2*(r_main/r_mag)^2*(5*z^2/r_mag^2-1);
    %               -mu*y/r_mag^3*J2*3/2*(r_main/r_mag)^2*(5*z^2/r_mag^2-1);
    %               -mu*z/r_mag^3*J2*3/2*(r_main/r_mag)^2*(3-5*z^2/r_mag^2)];
    
    J2_perturb = 3*J2*mu*r_main^2/(2*r_mag^5) * [(5*z^2/r_mag^2-1)*x;
                                                 (5*z^2/r_mag^2-1)*y;
                                                 (5*z^2/r_mag^2-3)*z];

    % EOM
    rddot = -mu*r/r_mag^3 + J2_perturb;
    zdot = [rdot; rddot];
end