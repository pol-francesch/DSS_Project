function zdot = meanDfq(t,z,rE,mu,J2)
    sma = z(1); ecc = z(2);
    ex  = z(3); ey  = z(4); inc = z(5);

    p_bar = sma*(1-ecc^2);
    n_bar = sqrt(mu/(sma^3*(1-ecc^2)^3));
    dRaan = -3/2*J2*(rE/p_bar)^2*n_bar*cos(inc);
    % dAop = 3/4*J2*(rE/p_bar)^2*n_bar*(5*cos(inc)^2-1);
    % eta_bar = sqrt(1-ecc^2);
    % dM0 = n_bar + 3/4*J2*(rE/p_bar)^2*n_bar*eta_bar*(3*cos(inc)^2-1);
    
    e_bar = 1-(ex^2+ey^2);
    du = 3/4*n_bar*J2*(rE/(sma*e_bar))^2*(sqrt(e_bar)*(3*cos(inc)^2-1) + (5*cos(inc)^2 - 1));
    dex = -3/4*n_bar*J2*(rE/(sma*e_bar))^2*ey*(5*cos(inc)^2 - 1);
    dey = 3/4*n_bar*J2*(rE/(sma*e_bar))^2*ex*(5*cos(inc)^2 - 1);
    
    % zdot = [0; 0; 0; dRaan; dAop; dM0];
    zdot = [0; 0; dex; dey; 0; dRaan; du];
end