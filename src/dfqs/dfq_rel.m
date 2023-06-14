function stateDot = dfq_rel(t,state,mu)
    % Calculates the derivative of z at t
    % Inputs:
    %   t        - current time, s
    %   state    - state vector [rho, rhodot, r0, theta0, rdot0, thetadot0]
    %   mu       - gravitational parameter of central body, m^3/s^2
    % Outputs:
    %   stateDot - derivative of state vector 
    %             [rhodot, rhoddot, rdot0, thetadot0, rddot0, thetaddot0]
    % ---------------------------------------------------------------------

    % Current state
    x = state(1); y = state(2); z = state(3);
    xdot = state(4); ydot = state(5); %zdot = state(6); (does not appear in EOM)
    rhodot = state(4:6); 
    r0 = state(7); % scalars from here down
    %theta0 = state(8); (does not appear in EOM)
    rdot0 = state(9);
    thetadot0 = state(10);

    % EOM of the chief
    rddot0 = r0*thetadot0^2 - mu/r0^2;
    thetaddot0 = -2*rdot0*thetadot0/r0;

    % Relative EOM of the deputy w.r.t. chief (in RTN frame)
    xddot = 2*thetadot0*ydot + thetaddot0*y + thetadot0^2*x - mu*(r0 + x)/((r0+x)^2+y^2+z^2)^(3/2) + mu/r0^2;
    yddot = -2*thetadot0*xdot - thetaddot0*x + thetadot0^2*y - mu*y / ((r0+x)^2+y^2+z^2)^(3/2);
    zddot = -mu*z/((r0+x)^2+y^2+z^2)^(3/2);
    rhoddot = [xddot; yddot; zddot];

    stateDot = [rhodot; rhoddot; rdot0; thetadot0; rddot0; thetaddot0];
end

