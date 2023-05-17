function qns_oe_d = roe2qns(qns_ref,roe)
    % Given the quasi-nonsingular orbital elements of the chief as reference and
    % non-dimensionalized relative orbit element offsets, produces the deputy QNS orbital elements
    % Inputs: qns_oe_c
    %         roe_init = [dsma;dlambda;dex;dey;dix;diy]
    % Output: [a, ex, ey, ix, iy, lambda]

    % Unpack ROEs
    dsma = roe(1); dlambda = roe(2); dex = roe(3);
    dey = roe(4);  dix = roe(5);     diy = roe(6);

    % Unpack OEs
    sma_c = qns_ref(1); ex_c = qns_ref(2);   ey_c = qns_ref(3); 
    inc_c = qns_ref(4); raan_c = qns_ref(5); u_c = qns_ref(6);
    
    % Mean argument of latitude: uM = v + w
    uM_c = true2mean(u_c, sqrt(ex_c^2 + ey_c^2));

    sma_d = sma_c + dsma*sma_c;
    ex_d = ex_c + dex;
    ey_d = ey_c + dey;
    inc_d = inc_c + dix;
    raan_d = raan_c + diy / sin(inc_c);
    uM_d = uM_c + dlambda - (raan_d - raan_c) * cos(inc_c);
    
    % True argument of latitude
    u_d = mean2true(uM_d, sqrt(ex_d^2 + ey_d^2), 1e-8);

    qns_oe_d = [sma_d, ex_d, ey_d, inc_d, raan_d, u_d];
end