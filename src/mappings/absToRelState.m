function [rho_RTN, rhodot_RTN] = absToRelState(mu,chief,deputy)
    % Calculates relative RTN position/velocity from absolute chief and
    % deputy positions in ECI
    % 
    % Inputs:
    %   mu          - gravitational parameter of central body [m^3/s^2]
    %   t           - time
    %   chief       - chief position and velocity [m, m/s]
    %   deputy      - deputy position and velocity [m, m/s]
    % Outputs:

    % Orbital elements
    oe_c = eci2oe(mu, chief)';
    qns_oe_c = singular2qns(oe_c);

    sma_c = oe_c(1); ecc_c = oe_c(2); inc_c = oe_c(3);
    raan_c = oe_c(4); u_c = qns_oe_c(6);

    % Keplerian
    % In ECI
    rho_ECI = (deputy(1:3)-chief(1:3));
    rhodot_ECI = (deputy(4:6)-chief(4:6));

    % Convert from ECI to RTN
    T_ECI2RTN = (rtn2eci_rot(inc_c,raan_c,u_c))'; % rotation matrix from ECI to RTN
    u_dot = sqrt(mu/(sma_c^3*(1-ecc_c^2)^3)) * (1+ecc_c*cos(u_c))^2;
    w_rtn2eci = [0;0;u_dot]; % rotation of the RTN frame with respect to ECI

    rho_RTN = T_ECI2RTN*rho_ECI;
    rhodot_RTN = T_ECI2RTN*rhodot_ECI - cross(w_rtn2eci,rho_RTN);

    rho_RTN = rho_RTN';
    rhodot_RTN = rhodot_RTN';
end

