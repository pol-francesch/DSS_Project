function f = ecc2true(E, e)
% ecc2true Computes true anomaly given eccentric anomaly and eccentricity
%
% Inputs:
%       E - eccentric anomaly [rad]
%       e - eccentricity of orbit
%
% Outputs:
%       f - mean anomaly [rad]

% E = mod(E, 2*pi);

f = atan2(sin(E)*sqrt(1-e^2),cos(E)-e);
% f = mod(atan2(sin(E)*sqrt(1-e^2),cos(E)-e),2*pi);
% f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));

end