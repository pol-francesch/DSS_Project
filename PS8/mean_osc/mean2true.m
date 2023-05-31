function f = mean2true(M, e, tol)
% mean2true solves Kepler's equation for true anomaly
%
%   Note: This function uses a Newton-Raphson method to numerically compute
%   the correct value for f
%
% Inputs:
%     M - mean anomaly [rad]
%     e - eccentricity of orbit
%   tol - tolerance for Newton-Raphson iterator
%
% Outputs:
%     f - true anomaly [rad]

if nargin == 2
    tol = 10e-10;
end

E = mean2ecc(M,e,tol);
f = ecc2true(E,e);

end