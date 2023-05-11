%% Reconfiguration simulation from mode 1 to mode 2 using least squares
% AA 279D Problem Set 5
% Sydney Hsu and Pol Francesch

clear; clc; close all;

% Constants
mu = 3.986e14;          % Earth gravitational parameter [m^3/s^2]
rE = 6378.127e3;        % Earth radius [m]
J2 = 0.00108263;        % Earth's J2 coefficient
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesFontSize',18)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 Part 6
% Initial Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial chief singular orbital parameters
sma_c  = 6892.927e3;      % semi-major axis [m]
ecc_c  = 1e-4;            % eccentricity component in I
inc_c  = deg2rad(97.44);  % inclination [rad]
raan_c = deg2rad(270);    % RAAN [rad]
aop_c  = deg2rad(0);      % aop [rad]
ta_c   = deg2rad(0);      % true anomaly [rad]

% Eccentricity vector
ex_c = ecc_c*cos(aop_c);
ey_c = ecc_c*sin(aop_c);
n = sqrt(mu/sma_c^3); % mean motion

oe_init_c = [sma_c, ecc_c, inc_c, raan_c, aop_c, ta_c]; % combine the above into a single vector
qns_oe_init_c = singular2qns(oe_init_c)'; % chief QNS OE (non-dimensionalized)

% Mode 1: Preliminary DEM generation offsets (phase C1)
dsma_1 = 0;
dlambda_1 = 0;
de_1 = 260 / sma_c;
di_1 = 222 / sma_c;
phase_diff_1 = 200; % designed phase difference % TN plane not relevant as it is safe in the other planes; along-track is the greatest uncertainty
theta_1 = deg2rad(90); % phase angle of inclination

% Calculate deputy initial orbit parameters
phi_1 = theta_1 - deg2rad(phase_diff_1); % phase angle of eccentricity
dex_1 = de_1*cos(phi_1);  dey_1 = de_1*sin(phi_1);
dix_1 = di_1*cos(theta_1); diy_1 = di_1*sin(theta_1);
% Initial conditions - quasi-nonsingular ROEs
roe_init1 = [dsma_1;dlambda_1;dex_1;dey_1;dix_1;diy_1]; % non-dimensionalized

% Initial pos/vel of deputy
qns_oe_init_d1 = roe2qns(qns_oe_init_c,roe_init1); % QNS orbital elements of the deputy
oe_init_d1 = qns2singular(qns_oe_init_d1); % singular orbital elements of the deputy [sma, ecc, inc, raan, aop, ta]
oe_init_d1_j2 = mean2osc(oe_init_d1, 1);
[r_init_d,v_init_d] = singular_oe2rv(mu,oe_init_d1_j2); % deputy position and velocity in inertial frame, J2 perturbed

% Initial pos/vel of chief
oe_init_c = qns2singular(qns_oe_init_c); % singular oe of chief
oe_init_c_j2 = mean2osc(oe_init_c, 1);
[r_init_c,v_init_c] = singular_oe2rv(mu,oe_init_c_j2); % chief position and velocity in inertial frame, J2-perturbed

% Orbital period of chief
T = 2*pi*sqrt(sma_c^3/mu); % [sec]


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 Part a
% Formation reconfig. parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define desired end state after config.
% Mode 2: phase C2
dsma_2 = 0;
dlambda_2 = 0;
de_2 = 297 / sma_c;
di_2 = 260 / sma_c;
phase_diff_2 = 210; % designed phase difference % TN plane not relevant as it is safe in the other planes; along-track is the greatest uncertainty
theta_2 = theta_1; % phase angle of inclination
% Calculate deputy initial orbit parameters
phi_2 = theta_2 - deg2rad(phase_diff_2); % phase angle of eccentricity
dex_2 = de_2*cos(phi_2);  dey_2 = de_2*sin(phi_2);
dix_2 = di_2*cos(theta_2); diy_2 = di_2*sin(theta_2);
% Initial conditions - quasi-nonsingular ROEs
roe_init2 = [dsma_2;dlambda_2;dex_2;dey_2;dix_2;diy_2];
droe1to2 = roe_init2 - roe_init1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 Part b
% Least squares solution to Reconfigure from 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 4; % Choose number of maneuvers
u_ks = [0, pi, 0, pi]; % selected maneuver locations
%delta_us = [0, pi, pi, pi]; % location differences (mean arg. of lat) between maneuvers
delta_us = [0, pi, 2*pi, 3*pi]+2*pi;
t_man = [0, T/2, T, 3/2*T]+T; % time of maneuvers

% Apply least-squares maneuver
[dV_lsqr, dV_lsqr_err] = least_squares_maneuvers(mu, J2, rE, N,u_ks,delta_us, droe1to2, oe_init_c);
dV = reshape(dV_lsqr,[3,N]);

% Verification with numerical propagation (FODE)
orb_revs = 5;
% days_elapsed = orb_revs*T/3600;
stepSize = T/1000;
tspan = 0:stepSize:T*orb_revs;
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances
rv_init_c = [r_init_c;v_init_c];
rv_init_d = [r_init_d;v_init_d];

%%
% In-plane
% dv_t = n/2*sma_c*dsma_2;
% dv_r = n*sma_c*sqrt(de_2^2+dsma_2^2);
% u_M_ip = atan2(dv_r,(2*dv_t)) - atan2(dey_2,dex_2);

de_1to2 = sqrt(droe1to2(3)^2+droe1to2(4)^2);
dv_t = n/2*sma_c*droe1to2(1);
dv_r = n*sma_c*sqrt(de_1to2^2+droe1to2(1)^2);
u_M_ip = atan2(dv_r,(2*dv_t)) - atan2(droe1to2(4),droe1to2(3));

% Out-of-plane
% dv_n = n*sma_c*di_2;
% u_M_oop = atan2(diy_2,dix_2);
di_1to2 = sqrt(droe1to2(5)^2+droe1to2(6)^2);
dv_n = n*sma_c*di_1to2;
u_M_oop = atan2(droe1to2(6),droe1to2(5));


%% 
% Propagate the chief and deputy orbits
[t,state_eci_c] = ode45(@(t,z) dfq_ptb(t,z,rE,mu,J2,0,0,0),tspan,rv_init_c,options);
% [~,state_eci_d] = ode45(@(t,z) dfq_ptb(t,z,rE,mu,J2,stepSize,t_man,dV),tspan,rv_init_d,options);
% [~,state_eci_d] = ode45(@(t,z) dfq_ptb(t,z,rE,mu,J2,0,0,0),tspan,rv_init_d,options);

%earth_plot(rE,state_eci_c,state_eci_d);

% Manually set dV
% dV = [[0, 0, n*sma_c*di_2]]';
% t_man = 3*T/2;

% Initialize state arrays
state_eci_d = []; % zeros(steps,6);

% Manually propagate to each maneuver time
t0 = 0;
r0_d = rv_init_d;
true_tspan = [];
man_dVs = dV;

for j=1:length(t_man)
    % Propagate to maneuver
    tspan_to_man = t0:stepSize:t_man(j);
    [~, deputy] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_to_man,r0_d,options);

    % Perform maneuver
    post_man = deputy(end,:) + [0,0,0,man_dVs(:,j)'];

    % Add to common state
    if tspan_to_man(end) == t_man(j)
        deputy = deputy(1:end-1,:);
        tspan_to_man = tspan_to_man(1:end-1);
    end
    % state_eci_d = [state_eci_d; deputy; post_man];
    % true_tspan  = [true_tspan, tspan_to_man, t_man(j)];
    state_eci_d = [state_eci_d; deputy];
    true_tspan  = [true_tspan, tspan_to_man];

    r0_d = post_man;
    t0   = t_man(j);
end

% Finish the propagation 
if t_man(end) < tspan(end)
    tspan_to_end = t0:stepSize:tspan(end);
    [~, deputy] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_to_end,r0_d,options);

    state_eci_d = [state_eci_d; deputy];
    true_tspan  = [true_tspan, tspan_to_end];
end

% tspan = true_tspan;
state_eci_d = interp1(true_tspan, state_eci_d, tspan);

earth_plot(rE,state_eci_c,state_eci_d);

%% Convert position and velocity to OE & ROE
steps = length(tspan);
% Osculating orbital elements
osc_qns_c = zeros(steps,6);
osc_sing_c = zeros(steps,6);
osc_qns_d = zeros(steps,6);
osc_qns_roe = zeros(steps,6);

osc_qns_j2 = zeros(steps, 6, 2);
mean_qns_j2 = zeros(steps, 6, 2);

osc_qns_j2_roe = zeros(steps, 6);
mean_qns_j2_roe = zeros(steps, 6);

for j=1:steps
    oe_c_sing = rv2singular_oe(mu,state_eci_c(j,1:3),state_eci_c(j,4:6)); % pos/vel to singular oe
    osc_sing_c(j,:) = oe_c_sing;
    qns_c = singular2qns(oe_c_sing); % singular oe to QNS oe;
    osc_qns_c(j,:) = qns_c;
    oe_d_sing = rv2singular_oe(mu,state_eci_d(j,1:3),state_eci_d(j,4:6)); % pos/vel to singular oe
    qns_d = singular2qns(oe_d_sing); % singular oe to QNS oe
    osc_qns_d(j,:) = qns_d;
    osc_qns_roe(j,:) = oes2roe(qns_c,qns_d);

    [osc_oe, mean_oe, osc_roe, mean_roe] = rv2osc_and_mean(mu, rE, J2, 1, 0, state_eci_c(j,:), state_eci_d(j,:));
    osc_qns_j2(j,:,:) = osc_oe';
    mean_qns_j2(j,:,:) = mean_oe';
    osc_qns_j2_roe(j,:) = osc_roe';
    mean_qns_j2_roe(j,:) = mean_roe';
end


% Plotting orbital elements
orbit_periods = t/T;
qns_labels = ["$a$ (m)", "$e_x$", "$e_y$", "$i$ (rad)",...
                "$\Omega$ (rad)", "$u$ (rad)"];
qns_roe_labels = ["$\delta a$ (m)", "$\delta \lambda$", "$\delta e_x$", ...
               "$\delta e_y$", "$\delta i_x$", "$\delta i_y$"];

figure('Name','Chief QNS Absolute')
plot_oe(orbit_periods,osc_qns_c,qns_labels)
figure('Name','Deputy QNS Absolute')
plot_oe(orbit_periods,osc_qns_d,qns_labels)
figure('Name','QNS ROE')
plot_oe(orbit_periods,osc_qns_roe,qns_roe_labels);
figure('Name','ROE State Space')
osc_qns_roe_norm = sma_c*osc_qns_roe;
plot_roe_state(osc_qns_roe_norm);
figure('Name','ROE State Space w/ init & final')
plot_roe_state_w_init_and_final(osc_qns_roe_norm, roe_init1*sma_c, roe_init2*sma_c);

% figure('Name','Singular OE Chief')
% plot_oe(orbit_periods,osc_sing_c,qns_labels);
%%
% figure('Name','Initial State Space')
% osc_qns_roe_norm = sma_c*osc_qns_roe;
% plot_roe_state(osc_qns_roe_norm(1:1340,:));
% 
% figure('Name','Final State Space')
% osc_qns_roe_norm = sma_c*osc_qns_roe;
% plot_roe_state(osc_qns_roe_norm(1340:end,:));


%% Perturbed absolute chief
% plot_osc_mean(orbit_periods, mean_qns_j2(:,:,1), osc_qns_j2(:,:,1), qns_labels);
% 
% % Perturbed absolute deputy
% plot_osc_mean(orbit_periods, mean_qns_j2(:,:,2), osc_qns_j2(:,:,2), qns_labels);
% 
% % Perturbed relative
% plot_osc_mean(orbit_periods, mean_qns_j2_roe, osc_qns_j2_roe, qns_roe_labels);




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 Part c
% Analytical single-impulse relative orbit control (D'Amico thesis, pg. 43)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % In-plane
% dv_t = n/2*sma_c*dsma_2;
% dv_r = n*sma_c*sqrt(de_2^2+dsma_2^2);
% u_M_ip = atan(dv_r/(2*dv_t)) - atan(dey_2/dex_2);
% 
% % Out-of-plane
% dv_n = n*sma_c*di_2;
% u_M_oop = atan(diy_2/dix_2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 
% Helper fcns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_roe_state_w_init_and_final(roe_set, roe_init, roe_final)
    % Relative eccentricity
    subplot(3,1,1); hold on; grid on; axis equal;
    plot(roe_set(:,3),roe_set(:,4),'.','MarkerSize',20);
    plot(roe_init(3), roe_init(4), '.', 'MarkerSize', 20);
    plot(roe_final(3), roe_final(4), '.', 'MarkerSize', 20);
    hold on
    xlabel('$a \delta e_x $ (m)')
    ylabel('$a \delta e_y$ (m)')

    % Relative inclination
    subplot(3,1,2); hold on; grid on; axis equal;
    plot(roe_set(:,5),roe_set(:,6),'.','MarkerSize',20); hold on
    plot(roe_init(5), roe_init(6), '.', 'MarkerSize', 20);
    plot(roe_final(5), roe_final(6), '.', 'MarkerSize', 20);
    hold on
    xlabel('$a \delta i_x $ (m)')
    ylabel('$a \delta i_y$ (m)')

    % Relative mean arg. of latitude and semi-major axis
    subplot(3,1,3); hold on; grid on; axis equal;
    plot(roe_set(:,2),roe_set(:,1),'.','MarkerSize',20); hold on
    plot(roe_init(2), roe_init(1), '.', 'MarkerSize', 20);
    plot(roe_final(2), roe_final(1), '.', 'MarkerSize', 20);
    hold on
    xlabel('$a \delta \lambda $ (m)')
    ylabel('$a \delta a$ (m)')
    legend('Time history', 'Initial', 'Final')
end
