function roe = singular_oe2roe(oe_c, oe_d)
    % Convert singular orbital elements to QNS ROEs 
    % Input: oe_c = [sma, e, inc, aop, raan, M]
    % Output: [sma, dlambda, dex, dey, dix, diy]
    sma_c = oe_c(1); e_c = oe_c(2); inc_c = oe_c(3); raan_c = oe_c(4); aop_c = oe_c(5); M_c = oe_c(6);
    sma_d = oe_d(1); e_d = oe_d(2); inc_d = oe_d(3); raan_d = oe_d(4); aop_d = oe_c(5); M_d = oe_d(6);

    roe = [ (sma_d-sma_c)/sma_c;
            (M_d + aop_d) - (M_c + aop_c) + (raan_d-raan_c)*cos(inc_c);
            e_d*cos(aop_d) - e_c*cos(aop_c);
            e_d*sin(aop_d) - e_c*sin(aop_c);
            inc_d-inc_c;
            (raan_d-raan_c)*sin(inc_c)]';

    % If an angle is close to 2pi, make it 0
    % for i=[2,5,6]
    %     if 2*pi - abs(roe(i)) < 1e-07
    %         roe(i) = 0;
    %     end
    % end

%     roe  = [(oe_d(1) - oe_c(1)) / oe_c(1);
%             (oe_d(6) - oe_c(6)) +  (oe_d(5) - oe_c(5))*cos(oe_c(4));
%             oe_d(2) - oe_c(2);
%             oe_d(3) - oe_c(3);
%             oe_d(4) - oe_c(4);
%             (oe_d(5) - oe_c(5))*sin(oe_c(4))]';

    % Make sure all angles are positive
    % roe(2) = wrapTo2Pi(roe(2));
    % roe(5:6) = wrapTo2Pi(roe(5:6));
end
