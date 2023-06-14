function [osc_sing_oe, mean_sing_oe, osc_qns_oe, mean_qns_oe, osc_roe, mean_roe] = rv2osc_and_mean(mu, rE, J2, J2_flag, deg_flag, chief, deputy)
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
    %   osc_qns_oe              - osculating quasi non-singular orbital
    %                         elements of chief and deputy
    %   mean_qns_oe             - osculating quasi non-singular orbital
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
    osc_qns_oe_chief = singular2qns(osc_sing_oe_chief); % [sma, ex, ey, inc, raan, u]
    mean_qns_oe_chief = singular2qns(mean_sing_oe_chief);

    % Deputy
    osc_qns_oe_deputy = singular2qns(osc_sing_oe_deputy);
    mean_qns_oe_deputy = singular2qns(mean_sing_oe_deputy);

    % QNS ROEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    osc_roe = oes2roe(osc_qns_oe_chief, osc_qns_oe_deputy);
    mean_roe = oes2roe(mean_qns_oe_chief, mean_qns_oe_deputy);
    % osc_roe = singular_oe2roe(osc_sing_oe_chief, osc_sing_oe_deputy);
    % mean_roe = singular_oe2roe(mean_sing_oe_chief, mean_sing_oe_deputy);

    % Degrees %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if deg_flag == 1
        disp("compute_osc_mean_oe: Degrees not yet implemented")
    end

    % Wrapping up %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Enforce row vectors as output
    osc_sing_oe = [reshape( osc_sing_oe_chief, [],1 ), reshape( osc_sing_oe_deputy, [],1 )];
    mean_sing_oe = [reshape( mean_sing_oe_chief, [],1 ), reshape( mean_sing_oe_chief, [],1 )];
    osc_qns_oe = [reshape( osc_qns_oe_chief, [],1 ), reshape( osc_qns_oe_deputy, [],1 )];
    mean_qns_oe = [reshape( mean_qns_oe_chief, [],1 ), reshape( mean_qns_oe_deputy, [],1 )];
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

    osc  = eci2oe(mu, rv)';
    mean = osc2mean(osc, J2_flag);
    % mean = osc;

    % Make sure all the angles are positive
    %mean(3:6) = wrapTo2Pi(mean(3:6));
    % mean(3:5) = wrapTo2Pi(mean(3:5));
end