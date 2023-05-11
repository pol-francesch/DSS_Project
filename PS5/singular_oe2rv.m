function [r_eci, v_eci] = singular_oe2rv (mu, oe)
    % OE2ECI Converts orbital elements to r , v in ECI frame
    %
    % Notes :
    % 1) In cases of equatorial and / or circular orbits , it is assumed
    % that valid orbital elements are provided as inputs ( ie . there is
    % no back - end validation )
    %
    % Inputs :
    %   mu - grav parameter [DU^3 / TU^2]
    %   oe - array of
    %       a - semi - major axis of orbit [DU]
    %       e - eccentricity of orbit
    %       i - inclination of orbit [rad]
    %       Om - right ascension of the ascending node [rad]
    %       w - argument of periapsis [rad]
    %       v - true anomaly [rad]
    %
    % Outputs :
    %   r_eci - 3 x1 vector of radius in ECI frame [DU]
    %   v_eci - 3 x1 vector of velocity in ECI frame [DU / TU]
    
    % Unpack
    a = oe(1); e = oe(2); i = oe(3);
    Om = oe(4); w = oe(5); v = oe(6);

    n = sqrt (mu / a^3) ; % rad / TU
    E = nu2E (v, e) ; % rad

    % Compute radius and velocity of orbit in perifocal coordinates
    r = a*(1-e^2)/(1+e*cos(v));
    rPeri = [r*cos(v); r*sin(v); 0];
    vPeri = sqrt(mu/(a*(1-e^2)))*[-sin(v); e+cos(v); 0];

    % Develop rotation matrix depending on orbit shape / inclination
    rotz_rad = @(x) rotz(rad2deg(x));
    rotx_rad = @(x) rotx(rad2deg(x));
    
    if i == 0 && e ~= 0 % Equatorial + elliptical
        rotPeri2ECI = rotz(w) ;
    elseif e == 0 && i ~= 0 % Circular + inclined
        rotPeri2ECI = rotz_rad(Om) * rotx_rad(i) ;
    elseif i == 0 && e == 0 % Equatorial + circular
        rotPeri2ECI = 1;
    else % Elliptical + inclined
        rotPeri2ECI = rotz_rad(Om) * rotx_rad(i) * rotz_rad(w) ;
    end

    % Rotate vectors into ECI frame
    r_eci = rotPeri2ECI * rPeri;
    v_eci = rotPeri2ECI * vPeri;
end

