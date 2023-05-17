function [osc_sing_c, osc_sing_d, osc_qns_c, osc_qns_d, mean_qns_c, mean_qns_d, osc_qns_roe, mean_qns_roe] = rv2oe_and_roe(mu, rE, J2, J2_flag, deg_flag, chief, deputy)
    steps = length(chief);

    osc_qns = zeros(steps, 6, 2);
    mean_qns = zeros(steps, 6, 2);
    osc_qns_roe = zeros(steps, 6);
    mean_qns_roe = zeros(steps, 6);
    sing_oe = zeros(steps,6,2);

    for j=1:steps
        [osc_sing_oe, osc_oe, mean_oe, osc_roe, mean_roe] = rv2osc_and_mean(mu, rE, J2, J2_flag, deg_flag, chief(j,:), deputy(j,:));
        osc_qns(j,:,:) = osc_oe';
        mean_qns(j,:,:) = mean_oe';
        osc_qns_roe(j,:) = osc_roe';
        mean_qns_roe(j,:) = mean_roe';
        sing_oe(j,:,:) = osc_sing_oe';
    end

    osc_qns_c = osc_qns(:,:,1);
    osc_qns_d = osc_qns(:,:,2);

    mean_qns_c = mean_qns(:,:,1);
    mean_qns_d = mean_qns(:,:,2);

    osc_sing_c = sing_oe(:,:,1);
    osc_sing_d = sing_oe(:,:,1);
end

