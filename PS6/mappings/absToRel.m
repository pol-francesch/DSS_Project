function [rho, rhodot] = absToRel(mu,t,chief,deputy) % not used: sma0,ecc0,inc0,raan0)
    % Calculates relative RTN position history from absolute chief and deputy positions
    % 
    % Inputs:
    %   mu          - gravitational parameter of central body [m^3/s^2]
    %   t           - time
    %   chief       - chief position and velocity [m, m/s]
    %   deputy      - deputy position and velocity [m, m/s]
    % Outputs:

    rtn = zeros(length(t),6);
    
    for j=1:length(t)
        % Orbital elements
        oe_c = eci2oe(mu, chief(j,:))';
        qns_oe_c = singular2qns(oe_c);

        sma_c = oe_c(1); ecc_c = oe_c(2); inc_c = oe_c(3);
        raan_c = oe_c(4); u_c = qns_oe_c(6);

        % Keplerian
        % In ECI
        rho_ECI = (deputy(j,1:3)-chief(j,1:3))';
        rhodot_ECI = (deputy(j,4:6)-chief(j,4:6))';
    
        % Convert from ECI to RTN
        T_ECI2RTN = (rtn2eci_rot(inc_c,raan_c,u_c))'; % rotation matrix from ECI to RTN
        u_dot = sqrt(mu/(sma_c^3*(1-ecc_c^2)^3)) * (1+ecc_c*cos(u_c))^2;
        w_rtn2eci = [0;0;u_dot]; % rotation of the RTN frame with respect to ECI
    
        rho_RTN = T_ECI2RTN*rho_ECI;
        rhodot_RTN = T_ECI2RTN*rhodot_ECI - cross(w_rtn2eci,rho_RTN);
        rtn(j,:) = [rho_RTN;rhodot_RTN]';
    end

    rho = rtn(:, 1:3);
    rhodot = rtn(:,4:6);
end