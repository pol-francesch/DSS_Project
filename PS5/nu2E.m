function E = nu2E ( v , e )
    % nu2E Computes eccentric anomaly given true anomaly and
    % eccentricity
    %
    % Inputs :
    % v - true anomaly [ rad ]
    % e - eccentricity of orbit
    %
    % Outputs :
    % E - eccentric anomaly [ rad ]

    %E = acos (( e + cos ( v ) ) /(1 + e * cos ( v ) ) ) ;
    E = 2*atan(sqrt((1-e)/(1+e))*tan(v/2));

    % Make sure E sits in the correct semi - plane
    if v > pi
        E = 2* pi - E ;
    end
end