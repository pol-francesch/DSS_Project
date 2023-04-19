function [r_0_RTN, rdot_0_RTN,thetadot0_0] = relative_motion(r,v,sma,ecc,inc,raan,u,mu) %r0,v0,r1,v1,sma,ecc,inc,raan,u,mu)
    % Takes in absolute position and velocity vectors
    % in inertial frame and converts to relative in the RTN frame 

    % Coriolis theorem
    thetadot0_0 = sqrt(mu/(sma^3*(1-ecc^2)^3)) * (1+ecc*cos(u))^2; % RTN coords / ECI frame
    w_RTNinECI = [0;0;thetadot0_0]; % angular velocity of RTN frame w.r.t. ECI expressed in the RTN frame
    T_ECI2RTN = (rtn2eci_rot(inc,raan,u))'; % rotation matrix from ECI to RTN
    r_0_RTN = T_ECI2RTN*r;
    rdot_0_RTN = T_ECI2RTN*v - cross(w_RTNinECI,r_0_RTN); % RTN coords / RTN frame

end