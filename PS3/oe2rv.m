function [r_ECI,v_ECI] = oe2rv(mu,sma,ex,ey,inc,raan,u)
    % Keplerian orbital elements to position and velocity in ECI
    % Inputs:
    %   mu          - gravitaional parameter [m^3/s^2]
    %   sma         - semi-major axis [m]
    %   ex          - eccentricity along I
    %   ey          - eccentricity along J
    %   inc         - inclination [rad]
    %   raan        - right ascencion of the ascending node [rad]
    %   u           - argument of latitude [rad]
    % Outputs:
    %   r_ECI       - position in ECI coordinate frame, m
    %   v_ECI       - velocity in ECI coordinate frame, m/s
     
    ecc = sqrt(ex^2 + ey^2);

    % Perifocal coordinates
    r_PQW = sma*(1-ecc^2)/(1+ecc*cos(u)) * [cos(u); sin(u); 0]; % position, m
    v_PQW = sqrt(mu/(sma*(1-ecc^2))) * [-sin(u); ecc+cos(u); 0]; % velocity, m/s
    
    % Perifocal to ECI rotation matrices
    R_x = @(a) [1,        0,        0;
                0,        cos(a),  sin(a);
                0,        -sin(a), cos(a)];
    R_z = @(a) [cos(a),  sin(a),  0;
                -sin(a), cos(a),  0;
                0,        0,        1];

    R_PQW_IJK = R_z(-raan)*R_x(-inc);
    
    % ECI coordinates
    r_ECI = R_PQW_IJK * r_PQW;
    v_ECI = R_PQW_IJK * v_PQW;
end

