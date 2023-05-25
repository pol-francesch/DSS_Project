function M = true2mean(f, e)
    % true2mean solves Kepler's equation for mean anomaly
    %
    % Inputs:
    %     f - true anomaly [rad]
    %     e - eccentricity of orbit
    %
    % Outputs:
    %     M - mean anomaly [rad]
    
    f = mod(f, 2*pi);
    
    E = 2*atan(sqrt((1-e)/(1+e)) * tan(f/2));
    
    % Make sure E sits in the correct semi-plane
    E = wrapTo2Pi(E);
    
    % OR
    % E = 2*atan(sqrt((1-e)/(1+e)) * tan(f/2));
    
    E = mod(E, 2*pi);
    
    % Kepler's equation
    M = E - e*sin(E);
end