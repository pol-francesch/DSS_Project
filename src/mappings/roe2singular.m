function oe_d = roe2singular(oe_c,roe)
    % Given the singular orbital elements of the chief as reference and
    % non-dimensionalized relative orbit element offsets, produces the deputy singular orbital elements
    % Inputs: oe_ref
    %         roe_init = [dsma;dlambda;dex;dey;dix;diy]
    % Output: [a, e, i, Omega, omega, M]

    % Unpack ROEs
    dsma = roe(1); dlambda = roe(2); dex = roe(3);
    dey = roe(4);  dix = roe(5);     diy = roe(6);

    % Unpack OEs
    sma_c = oe_c(1);  e_c = oe_c(2); inc_c = oe_c(3); 
    raan_c = oe_c(4); w_c = oe_c(5); M_c = oe_c(6);

    sma_d = sma_c + dsma*sma_c;
    inc_d = inc_c + dix;
    raan_d = raan_c + diy / sin(inc_c);
    
    w_d = atan2(dey + e_c*sin(w_c), dex + e_c*cos(w_c));
    M_d = dlambda - (raan_d - raan_c)*cos(inc_c) + (M_c + w_c) - w_d;
    e_d = (dex + e_c*cos(w_c) - dey - e_c*sin(w_c)) / (cos(w_d) - sin(w_d));

    oe_d = [sma_d, e_d, inc_d, raan_d, w_d, M_d];

end

