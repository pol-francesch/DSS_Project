function E = calc_eccentric(M, ecc, eps)
    % Calculates eccentric anomaly using Newton-Raphson algorithm
    % Inputs:
    %   M       - mean anomaly, rad
    %   ecc     - eccentricity
    %   eps     - tolerance, rad
    % Outputs:
    %   E       - eccentric anomaly, rad
    
    % Parameters for iteration
    E = pi;
    max_iter = 1000;
    iter = 0;
    delta = 1000;
    
    % Newton-Raphson algorithm
    while abs(delta) > eps && iter < max_iter
        delta = - (E - ecc*sin(E) - M) / (1 - ecc*cos(E));

        % Update parameters
        E = E + delta;
        iter = iter + 1;
    end
end

