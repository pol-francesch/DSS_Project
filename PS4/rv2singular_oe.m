function oe = rv2singular_oe(mu, r_ijk, v_ijk)
    %{
        ECI2OE converts position , r , and velocity , v , vectors in the ECI frame
        to
        orbital elements .
        Inputs :
            mu - gravitational parameter of central body [DU^3 / TU^2]
            r_ijk - 3 x1 vector of radius in ECI frame [DU]
            v_ijk - 3 x1 vector of velocity in ECI frame [DU / TU]
        Outputs :
            oe = array with the following fields :
                a - semi - major axis of orbit [DU]
                e - eccentricity of orbit
                i - inclination of orbit [rad]
                Om - right ascension of the ascending node [rad]
                w - argument of periapsis [rad]
                anom - true anomaly [rad]
    %}

    v = norm (v_ijk);
    r = norm (r_ijk);

    % Create angular momentum , mean motion , and eccentricity vectors :
    hVec = cross(r_ijk, v_ijk);
    h = norm(hVec);
    nVec = cross([0, 0, 1], hVec);
    n = norm(nVec);
    eVec = (1/mu) * ((v^2 - mu/r)*r_ijk - dot(r_ijk, v_ijk)*v_ijk);
    
    e = norm ( eVec );

    % Compute orbit shape / size :
    mechEnergy = 0.5*v^2 - mu/r;

    if e ~= 1 % Elliptical , Circular , Hyperbolic orbits
        a = - mu/(2*mechEnergy);
    else % Parabolic orbit
        a = inf;
    end

    % Compute the orientation of the orbit
    % Rounding to 12 decimal places to avoid imag numbers
    i = acos(round(hVec(3) / h, 12));
    Om = acos(round(nVec(1) / n, 12));
    w = acos(round(dot(nVec, eVec) / (n*e), 12));
    anom = acos(round(dot(eVec, r_ijk) / (e*r), 12));

    % Place angles in the correct range :
    if nVec(2) < 0
        Om = 2*pi - Om ;
    end
    if eVec(3) < 0
        w = 2*pi - w ;
    end
    if dot(r_ijk , v_ijk) < 0
        anom = 2*pi - anom ;
    end
    
    i = wrapTo2Pi(i); Om = wrapTo2Pi(Om);
    w = wrapTo2Pi(w); anom = wrapTo2Pi(anom);

    % If the angles are close to 2pi, make them 0
    if 2*pi - abs(i) < 1e-12
        i = 0;
    end
    if 2*pi - abs(Om) < 1e-12
        Om = 0;
    end
    if 2*pi - abs(w) < 1e-12
        w = 0;
    end
    if 2*pi - abs(anom) < 1e-12
        anom = 0;
    end

    % Assign everything to array
    oe = [a, e, i, Om, w, anom];
end

