function E = nu2E (v , e )
    % nu2E Computes eccentric anomaly given true anomaly and
    % eccentricity
    %
    % Inputs :
    % v - true anomaly [ rad ]
    % e - eccentricity of orbit
    %
    % Outputs :
    % E - eccentric anomaly [rad ]
    E = acos (( e + cos(v) ) /(1 + e * cos(v) ) ) ;

    % Make sure E sits in the correct semi - plane
    if v > pi
        E = 2* pi - E ;
    end
end

