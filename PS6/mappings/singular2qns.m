function qns_oe = singular2qns(oe)
    % Singular OE to QNS OE
    %
    % Inputs :
    %   oe - array of
    %       a - semi - major axis of orbit [DU]
    %       e - eccentricity of orbit
    %       i - inclination of orbit [rad]
    %       Om - right ascension of the ascending node [rad]
    %       w - argument of periapsis [rad]
    %       M - mean anomaly [rad]
    %
    % Outputs :
    %   qns_oe - array of
    %       a - semi - major axis of orbit [DU]
    %       ex - eccentricity of orbit along x
    %       ey - eccentricity of orbit along y
    %       i - inclination of orbit [rad]
    %       Om - right ascension of the ascending node [rad]
    %       u - mean argument of latitude [rad]

    % Unpack
    sma = oe(1); ecc = oe(2);  inc = oe(3);
    raan = oe(4);   aop = oe(5); M = oe(6);
    
    % Convert Singular OE to QNS OE
    qns_oe = [sma;
              ecc*cos(aop);
              ecc*sin(aop);
              inc;
              raan;
              wrapTo2Pi(M + aop)]';
end


