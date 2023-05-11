function zdot = dfq_ptb(t,z,r_main,mu,J2,tstep,t_man,dV_man)
    % Calculates the derivative of z at t with J2 (FODE)
    % Inputs:
    %   t       - current time, s
    %   z       - state vector [rx, ry, rz, vx, vy, vz] [m, m/s]
    %   r_main  - main body radius, m
    %   mu      - gravitaional parameter of central body, m^3/s^2
    %   J2      - J2 coefficient
    %   tstep   - time step of simulation
    %   t_man    - list of maneuver times
    %   dV_man - 3xN 
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

    % Total perturbation
    d = J2_perturb;

    % Check if we are within a time step of a maneuver if specified
    if tstep > 0
        for j = 1:length(t_man)
            t_diff = t_man(j) - t;
            if t_diff >= 0 && t_diff < 1 %tstep % if we are within a time step of a maneuver
                % Execute the maneuver

                % ECI to RTN
                % Coriolis theorem
                % thetadot0_0 = sqrt(mu/(sma^3*(1-ecc^2)^3)) * (1+ecc*cos(u))^2; % RTN coords / ECI frame
                % w_RTNinECI = [0;0;thetadot0_0]; % angular velocity of RTN frame w.r.t. ECI expressed in the RTN frame
                % T_ECI2RTN = (rtn2eci_rot(inc,raan,u))'; % rotation matrix from ECI to RTN
                % r_0_RTN = T_ECI2RTN*r;
                % rdot_0_RTN = T_ECI2RTN*v - cross(w_RTNinECI,r_0_RTN); % RTN coords / RTN frame

                d = J2_perturb + dV_man(:,j) / tstep; % add acceleration of the maneuver
            end
        end
    end
    
    % EOM
    rddot = -mu*r/r_mag^3 + d;
    zdot = [rdot; rddot];
end