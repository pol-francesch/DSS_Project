function [osc_oe,mean_oe,osc_roe,mean_roe] = rv2osc_and_mean(mu, rE, J2, J2_flag, deg_flag, chief, deputy)
    % Given the chief and deputy position and velocity we compute
    % the osculating and mean orbital elements for both, and the osculating
    % and mean relative orbital elements
    % 
    % Inputs:
    %   mu          - gravitational parameter of central body [m^3/s^2]
    %   rE          - earth radius [m]
    %   J2          - Earth J2 coefficient
    %   J2_flag     - whether J2 needs to be considered
    %   deg_flag    - whether results should be in degrees
    %   chief       - chief position and velocity [m, m/s]
    %   deputy      - deputy position and velocity [m, m/s]
    % Outputs:
    %   osc_oe              - osculating quasi non-singular orbital
    %                         elements of chief and deputy
    %   mean_oe             - osculating quasi non-singular orbital
    %                         elements of chief and deputy
    %   osc_roe             - osculating quasi non-singular relative 
    %                         orbital elements
    %   mean_roe            - osculating quasi non-singular relative 
    %                         orbital elements
    
    % Singular OEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Chief
    [osc_sing_oe_chief, mean_sing_oe_chief] = rv2singular_osc_and_mean(mu, rE, J2, J2_flag, chief);

    % Deputy
    [osc_sing_oe_deputy, mean_sing_oe_deputy] = rv2singular_osc_and_mean(mu, rE, J2, J2_flag, deputy);

    % QNS OEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Chief
    osc_qns_oe_chief = singular2qns(osc_sing_oe_chief);
    mean_qns_oe_chief = singular2qns(mean_sing_oe_chief);

    % Deputy
    osc_qns_oe_deputy = singular2qns(osc_sing_oe_deputy);
    mean_qns_oe_deputy = singular2qns(mean_sing_oe_deputy);

    % QNS ROEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    osc_roe = oes2roe(osc_qns_oe_chief, osc_qns_oe_deputy);
    mean_roe = oes2roe(mean_qns_oe_chief, mean_qns_oe_deputy);

    % Degrees %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if deg_flag == 1
        disp("compute_osc_mean_oe: Degrees not yet implemented")
    end

    % Wrapping up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    osc_oe = [osc_qns_oe_chief; osc_qns_oe_deputy];
    mean_oe = [mean_qns_oe_chief; mean_qns_oe_deputy];
end

function [osc, mean] = rv2singular_osc_and_mean(mu, rE, J2, J2_flag, rv)
    % Given the position and velocity we compute
    % the singular osculating and mean orbital elements 
    % 
    % Inputs:
    %   mu          - gravitational parameter of central body [m^3/s^2]
    %   rE          - earth radius [m]
    %   J2          - Earth J2 coefficient
    %   J2_flag     - whether J2 needs to be considered
    %   rv          - position and velocity [m, m/s]
    % Outputs:
    %   osc              - osculating singular orbital
    %                         elements of chief and deputy
    %   mean             - osculating singular orbital
    %                         elements of chief and deputy

    osc  = rv2singular_oe(mu, rv(1:3), rv(4:6));
    mean = osc2mean(osc, J2_flag);

    % Make sure all the angles are positive
    mean(3:6) = wrapTo2Pi(mean(3:6));
end

function roe = oes2roe(oe_c, oe_d)
    % oe_c = [sma, ex, ey, inc, aop, u]
    % QNS OEs to ROEs [sma, dlambda, dex, dey, dix, diy]
    roe  = [(oe_d(1) - oe_c(1)) / oe_c(1);
            (oe_d(6) - oe_c(6)) + (oe_d(5) - oe_c(5))*cos(oe_c(4));
            oe_d(2) - oe_c(2);
            oe_d(3) - oe_c(3);
            oe_d(4) - oe_c(4);
            (oe_d(5) - oe_c(5))*sin(oe_c(4))]';

    % Make sure all angles are positive
    roe(2) = wrapTo2Pi(roe(2));
    roe(5:6) = wrapTo2Pi(roe(5:6));
end