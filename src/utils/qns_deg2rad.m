function qns_oe = qns_deg2rad(qns_oe)
    %QNS_DEG2RAD Takes QNS OE in radians and converts the angles to deg
    % Make sure it is (N, 6)
    qns_oe(:, 4:6) = rad2deg(qns_oe(:, 4:6));
end

