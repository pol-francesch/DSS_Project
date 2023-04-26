function M = ta2ma(aop,ecc)
    % True anomaly to mean anomaly
    % Inputs: arg. of periapsis, eccentricity
    E = nu2E(aop,ecc); % eccentric anomaly
    M = E-ecc*sin(E);
end