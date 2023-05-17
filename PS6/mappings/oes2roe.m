function roe = oes2roe(oe_c, oe_d)
    % Convert quasi-nonsingular orbital elements to QNS ROEs 
    % Input: oe_c = [sma, ex, ey, inc, raan, u]
    % Output: [sma, dlambda, dex, dey, dix, diy]
    sma_c = oe_c(1); ex_c = oe_c(2); ey_c = oe_c(3); inc_c = oe_c(4); raan_c = oe_c(5); uM_c = oe_c(6);
    sma_d = oe_d(1); ex_d = oe_d(2); ey_d = oe_d(3); inc_d = oe_d(4); raan_d = oe_d(5); uM_d = oe_d(6);

    % uM_c = true2mean(u_c, sqrt(ex_c^2 + ey_c^2));
    % uM_d = true2mean(u_d, sqrt(ex_d^2 + ey_d^2));

    roe = [ (sma_d-sma_c)/sma_c;
            % wrapTo2Pi((uM_d-uM_c)+(raan_d-raan_c)*cos(inc_c));
            (uM_d-uM_c)+(raan_d-raan_c)*cos(inc_c);
            ex_d-ex_c;
            ey_d-ey_c;
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