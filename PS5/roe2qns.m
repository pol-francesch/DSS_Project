function qns_oe_d = roe2qns(qns_ref,roe)
    % Given the quasi-nonsingular orbital elements of the chief as reference and
    % non-dimensionalized relative orbit element offsets, produces the deputy QNS orbital elements
    % Inputs: qns_oe_c
    %         roe_init = [dsma;dlambda;dex;dey;dix;diy]
    % Output: [a, ex, ey, ix, iy, lambda]

    % Unpack ROEs
    dsma = roe(1);
    dlambda = roe(2);
    dex = roe(3);
    dey = roe(4);
    dix = roe(5);
    diy = roe(6);

    qns_oe_d = qns_ref;
    qns_oe_d(1) = qns_ref(1) + qns_ref(1)*dsma/qns_ref(1);
    qns_oe_d(2) = qns_ref(2) + dex; %/qns_ref(1);
    qns_oe_d(3) = qns_ref(3) + dey; %/qns_ref(1);
    qns_oe_d(4) = qns_ref(4) + dix; %/qns_ref(1);
    qns_oe_d(5) = qns_ref(5) + diy/sin(qns_ref(4)); % diy/qns_ref(1)/sin(qns_ref(4));
    qns_oe_d(6) = qns_ref(6) + dlambda - (qns_oe_d(5) - qns_ref(5))*cos(qns_ref(4));
end