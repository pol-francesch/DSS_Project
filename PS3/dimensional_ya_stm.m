function [stm] = dimensional_ya_stm(t,mu,sma,ecc,ta)
    % Calculates dimensional STM from YA solution
    % Inputs:
    % Outputs:
    
    n = sqrt(mu/sma^3);

    % Auxiliary variables
    k = 1 + ecc*cos(ta);
    eta = sqrt(1 - ecc^2);
    tau = n*t / eta^3;
    kp = - ecc*sin(ta);
    
    % STM
    mat1 = [1/k + 3/2*kp*tau, sin(ta), cos(ta), 0, 0, 0;
           -3/2*k*tau, (1 + 1/k)*cos(ta), -(1 + 1/k)*sin(ta), 1/k, 0, 0;
           0, 0, 0, 0, 1/k*sin(ta), 1/k*cos(ta);
           kp/2 - 3/2*k^2*(k-1)*tau, k^2*cos(ta), -k^2*sin(ta), 0, 0, 0;
           -3/2*(k+k^2*kp*tau), -(k^2 + 1)*sin(ta), -ecc - (k^2+1)*cos(ta), -kp, 0, 0;
           0, 0, 0, 0, ecc + cos(ta), -sin(ta)];
    mat2 = [sma*eta^2*eye(3), zeros(3);
            zeros(3),         sma*n/eta*eye(3)];
    stm = mat2 * mat1;
end

