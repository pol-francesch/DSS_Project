function [osc, mean] = rv2singular_osc_and_mean(mu, rE, J2, J2_flag, rv)
    % Given the position and velocity we compute
    % the singular osculating and mean orbital elements 
    % 
    % Inputs:
    %   mu          - gravitational parameter of central body [m^3/s^2]
    %   rE          - earth radius [m]
    %   J2          - Earth J2 coefficient
    %   J2_flag     - whether J2 needs to be considered
    %   rv          - position and velocity [m, m/s] in inertial frame (ECI)
    % Outputs:
    %   osc              - osculating singular orbital
    %                         elements of chief and deputy
    %   mean             - osculating singular orbital
    %                         elements of chief and deputy

    osc  = eci2oe(mu, rv)';
    mean = osc2mean(osc, J2_flag);
    % mean = osc;

    % Make sure all the angles are positive
    %mean(3:6) = wrapTo2Pi(mean(3:6));
    % mean(3:5) = wrapTo2Pi(mean(3:5));
end