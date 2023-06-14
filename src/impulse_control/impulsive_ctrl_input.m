function B = impulsive_ctrl_input(n,eta,oe_sing)
    % Computes the control input matrix for impulsive control given the
    % singular ROE state
    % Inputs: n - mean motion
    %         eta - auxiliary variable
    %         oe_sing - singular orbital elements of the chief [sma, ecc, inc, raan, aop, Ma]
    % Outputs: B - control input matrix to be multiplied by control acc. in RTN (B*u_RTN)

    % Unpack singular orbital elements
    sma0 = oe_sing(1); ecc0 = oe_sing(2); inc0 = oe_sing(3); aop0 = oe_sing(5); M0 = oe_sing(6);

    % Compute true anomaly (from Steindorf pg. 22)
    tol = 1e-8;
    E = mean2ecc(M0,ecc0,tol);
    f0 = ecc2true(E,ecc0);
    % True argument of latitude
    u0 = f0 + aop0;
    
    % Quasi-nonsingular eccentricity components
    ex0 = ecc0*cos(aop0);
    ey0 = ecc0*sin(aop0);

    B = (1/(sma0*n)) * [2/eta*ecc0*sin(f0), 2/eta*(1+ecc0*cos(f0)), 0; ...
                        -2*eta^2/(1+ecc0*cos(f0)), 0, 0;
                         eta*sin(u0), eta*(2+ecc0*cos(f0)*cos(u0)+ex0)/(1+ecc0*cos(f0)),  eta*ey0/tan(inc0)*(sin(u0)/(1+ecc0*cos(f0)));
                        -eta*cos(u0), eta*(2+ecc0*cos(f0)*sin(u0)+ey0)/(1+ecc0*cos(f0)), -eta*ex0/tan(inc0)*(sin(u0)/(1+ecc0*cos(f0)));
                        0, 0, eta*cos(u0)/(1+ecc0*cos(f0));
                        0, 0, eta*sin(u0)/(1+ecc0*cos(f0)) ];

end