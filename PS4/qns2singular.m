function oe = qns2singular(qns_oe)
    % Singular OE to QNS OE
    %
    % Inputs :
    %   qns_oe - array of
    %       a - semi - major axis of orbit [DU]
    %       ex - eccentricity of orbit along x
    %       ey - eccentricity of orbit along y
    %       i - inclination of orbit [rad]
    %       Om - right ascension of the ascending node [rad]
    %       u - true argument of latitude [rad]
    % Outputs :
    %   oe - array of
    %       a - semi - major axis of orbit [DU]
    %       e - eccentricity of orbit
    %       i - inclination of orbit [rad]
    %       Om - right ascension of the ascending node [rad]
    %       w - argument of periapsis [rad]
    %       v - true anomaly [rad]

    % Unpack
    sma = qns_oe(1); ex = qns_oe(2);   ey = qns_oe(3);
    inc = qns_oe(4); raan = qns_oe(5); u = qns_oe(6);

    % Convert QNS OE to Singular OE
    oe = [sma;
          sqrt(ex^2 + ey^2);
          inc;
          raan;
          atan2(ey, ex);
          u - atan2(ey,ex)]';
end
