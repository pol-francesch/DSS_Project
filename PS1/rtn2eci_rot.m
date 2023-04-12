function rot = rtn2eci_rot(inc, raan, aop, ta)
    % Rotation matrix from RTN to ECI
    % Inputs:
    %   a       - semi-major axis, km
    %   e       - eccentricity
    %   i       - inclination, rad
    %   raan    - right ascencion of the ascending node, rad
    %   aop     - argument of periapsis, rad
    %   ta      - true anomaly, rad
    % Outputs:
    %   rot     - rotation matrix

    R_z = @(a) [cos(a),  sin(a),  0;
                -sin(a), cos(a),  0;
                0,       0,       1];
    R_x = @(a) [1,        0,        0;
                0,        cos(a),  sin(a);
                0,        -sin(a), cos(a)];

    rot = R_z(-raan)*R_x(-inc)*R_z(-aop)*R_z(-ta);
end

