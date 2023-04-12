function [r_ECI,v_ECI] = oe2rv(mu,sma,ecc,inc,raan,aop,ta)
    % Keplerian orbital elements to position and velocity in ECI
    % Inputs:
    %   mu          - gravitaional parameter, m^3/s^2
    %   sma         - semi-major axis, m
    %   ecc         - eccentricity
    %   inc         - inclination, rad
    %   raan        - right ascencion of the ascending node, rad
    %   aop         - argument of periapsis, rad
    %   ta          - true anomaly, rad
    % Outputs:
    %   r_ECI       - position in ECI coordinate frame, m
    %   v_ECI       - velocity in ECI coordinate frame, m/s
   
    % Mean motion
    n = sqrt(mu/sma^3);   % 1/s
    
    % Eccentric anomaly
    E = 2*atan(sqrt((1-ecc)/(1+ecc)) * tan(ta/2));    % rad
    
    % Perifocal coordinates
    r_PQW = sma*(1-ecc^2)/(1+ecc*cos(ta)) * [cos(ta); sin(ta); 0]; % position, m
    v_PQW = sqrt(mu/(sma*(1-ecc^2))) * [-sin(ta); ecc+cos(ta); 0]; % velocity, m/s
    
    % Perifocal to ECI rotation matrices
    R_x = @(a) [1,        0,        0;
                0,        cos(a),  sin(a);
                0,        -sin(a), cos(a)];
    R_z = @(a) [cos(a),  sin(a),  0;
                -sin(a), cos(a),  0;
                0,        0,        1];
    R_PQW_IJK = R_z(-raan)*R_x(-inc)*R_z(-aop);
    
    % ECI coordinates
    r_ECI = R_PQW_IJK * r_PQW;
    v_ECI = R_PQW_IJK * v_PQW;
end

