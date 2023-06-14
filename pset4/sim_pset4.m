%% These are Relative Orbits!
% AA 279D Problem Set 4
% Sydney Hsu and Pol Francesch
clear; clc; close all;

initial_conditions_pset4;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part c
% Osculating & Mean quasi-nonsingular OE & ROE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numerical integration 
% Settings
orb_rev = 15;
stepSize = T/1000;
tspan = 0:stepSize:T*orb_rev;
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances
[r,v] = singular_oe2rv(mu,oe_init_c);
z0_chief = [r;v];
[r,v] = singular_oe2rv(mu,oe_init_d1);
z0_deputy = [r;v];

z0_chief_j2 = rv_init_c;
z0_deputy_j2 = rv_init_d;

% Numerical integration
[t,z_c_unptb] = ode45(@(t,z) dfq(t,z,rE,mu,0),tspan,z0_chief,options);
[~,z_d_unptb] = ode45(@(t,z) dfq(t,z,rE,mu,0),tspan,z0_deputy,options);

[~,z_c_j2] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan,z0_chief_j2,options);
[~,z_d_j2] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan,z0_deputy_j2,options);

% To OE and ROE
[~, ~, osc_qns_unptb_c, osc_qns_unptb_d, mean_qns_unptb_c, mean_qns_unptb_d, osc_qns_unptb_roe, mean_qns_unptb_roe] = rv2oe_and_roe(mu, rE, J2, 0, 0, z_c_unptb, z_d_unptb);

[~, ~, osc_qns_j2_c, osc_qns_j2_d, mean_qns_j2_c, mean_qns_j2_d, osc_qns_j2_roe, mean_qns_j2_roe] = rv2oe_and_roe(mu, rE, J2, 1, 0, z_c_j2, z_d_j2);

% Plot
orbit_periods = t/T;
abs_ylabels = ["$a$ (m)", "$e_x$", "$e_y$", "$i$ (deg)",...
                "$\Omega$ (deg)", "$u$ (deg)"];
rel_ylabels = ["$a \delta a (m)$ ", "$ a \delta \lambda (m)$", "$a \delta e_x (m)$", ...
               "$a \delta e_y (m)$", "$a \delta i_x (m)$", "$a \delta i_y (m)$"];

% Unperturbed absolute chief
figure('Name', 'Unperturbed absolute chief');
plot_osc_mean(orbit_periods, qns_deg2rad(mean_qns_unptb_c), qns_deg2rad(osc_qns_unptb_c), abs_ylabels);

% Unperturbed absolute deputy
figure('Name', 'Unperturbed absolute deputy');
plot_osc_mean(orbit_periods, qns_deg2rad(mean_qns_unptb_d), qns_deg2rad(osc_qns_unptb_d), abs_ylabels);

% Unperturbed relative
figure('Name', 'Unperturbed relative');
plot_osc_mean(orbit_periods, sma_c*mean_qns_unptb_roe, sma_c*osc_qns_unptb_roe, rel_ylabels);

% Perturbed absolute chief
figure('Name', 'Perturbed absolute chief');
plot_osc_mean(orbit_periods, qns_deg2rad(mean_qns_j2_c), qns_deg2rad(osc_qns_j2_c), abs_ylabels);

% Perturbed absolute deputy
figure('Name', 'Perturbed absolute deputy');
plot_osc_mean(orbit_periods, qns_deg2rad(mean_qns_j2_d), qns_deg2rad(osc_qns_j2_d), abs_ylabels);

% Perturbed relative
figure('Name', 'Perturbed relative');
plot_osc_mean(orbit_periods, sma_c*mean_qns_j2_roe, sma_c*osc_qns_j2_roe, rel_ylabels);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part 4
% RTN frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get data
[rho2_unptb, ~] = absToRel(mu, t, z_c_unptb, z_d_unptb);
[rho2_j2, ~] = absToRel(mu, t, z_c_j2, z_d_j2);

% Plot
figure('Name', 'RTN Unperturbed and J2')
% rtn_plot_pos([rho2_unptb(:,1), rho2_j2(:,1)], ...
%              [rho2_unptb(:,2), rho2_j2(:,2)], ...
%              [rho2_unptb(:,3), rho2_j2(:,3)])
plot_rtn_pos(rho2_unptb(:,1:3), '-', 2); hold on;
plot_rtn_pos(rho2_j2(:,1:3), '-', 2);
legend('Unperturbed','J2')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part 5
% Relative orbital element state space plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Osculating ROEs normalized by semi-major axis
osc_qns_unptb_roe_state = sma_c*osc_qns_unptb_roe; % units: m
osc_qns_j2_roe_state = sma_c*osc_qns_j2_roe; % units: m

% Mean ROEs normalized by semi-major axis
mean_qns_unptb_roe_state = sma_c*mean_qns_unptb_roe;
mean_qns_j2_roe_state = sma_c*mean_qns_j2_roe;

figure('Name','Unperturbed ROE State')
plot_roe_state_2(osc_qns_unptb_roe_state,mean_qns_unptb_roe_state)
legend('Osculating','Mean')

figure('Name','J2 ROE State')
plot_roe_state_2(osc_qns_j2_roe_state,mean_qns_j2_roe_state)
legend('Osculating','Mean')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part 6, 7
% New Initial ROEs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_i = -30/sma_c; % desired inclination change
u_M = 0; % loc. of maneuver at ascending node
n = sqrt(mu/sma_c^3);
deltav_N = delta_i * n * sma_c / cos(u_M); % [m/s]

% New initial conditions for relative motion [dimensions are m -> sma_c*dx]
dsma_new = 0;
dlambda_new = 100;
dex_new = 50;
dey_new = 100;
dix_new = 0; %30;
diy_new = 200;

roe_init2 = [dsma_new;dlambda_new;dex_new;dey_new;dix_new;diy_new] ./ sma_c; % non-dimensionalized

% Initial pos/vel of deputy
oe_init_d1 = roe2singular(oe_init_c,roe_init2);
qns_oe_init_d1 = singular2qns(oe_init_d1);

oe_init_d1_j2 = mean2osc(oe_init_d1, 1);

% Position and velocity
[r_init_d,v_init_d] = singular_oe2rv(mu,oe_init_d1_j2); % deputy position and velocity in inertial frame, J2 perturbed
rv_init_d = [r_init_d;v_init_d];

% Numerical integration
[t,z_d_j2_nodrift] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan,rv_init_d,options);

% Convert to osculating oe
[~, ~, osc_qns_nodrift_c, osc_qns_nodrift_d, mean_qns_nodrift_c, mean_qns_nodrift_d, osc_qns_nodrift_roe, mean_qns_nodrift_roe] = rv2oe_and_roe(mu, rE, J2, 1, 0, z_c_j2, z_d_j2_nodrift);

% Visualization of relative drift
figure('Name','Drift Correction ROE')
orbit_periods = t/T;
for i = 1:1:6
    subplot(2,3,i); hold on;
    plot(orbit_periods, osc_qns_nodrift_roe(:,i,1));
    xlabel('Orbital periods')
    ylabel(rel_ylabels(i));
    ytickformat('%.3f')
    grid on; hold off;
end

% State space plots
osc_state_roe = sma_c*osc_qns_nodrift_roe;
plot_roe_state_1(osc_state_roe)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part 8
% J2 STM for ROE and comparison to numerical integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions - quasi-nonsingular ROEs
state0 = roe_init1*sma_c; % units: [m]
state1_0 = [dsma_new;dlambda_new;dex_new;dey_new;dix_new;diy_new];

% Set up propagation
% orb_rev = 1000;
stepSize = T/1000;
tspan = 0:stepSize:T*orb_rev;
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances
state_j2 = zeros(length(tspan),6);
state_j2_nodrift = zeros(length(tspan),6);
for j=1:length(tspan)
    t = tspan(j);
    % Relative position and velocity
    state_j2(j,:) = stm_qns_roe_j2(t,sma_c,ecc_c,inc_c,aop_c)*state0;
    state_j2_nodrift(j,:) = stm_qns_roe_j2(t, sma_c, ecc_c, inc_c, aop_c)*state1_0;
end

figure('Name','J2 STM ROE State')
plot_roe_state_1(state_j2)

figure('Name', 'J2 STM ROE State New ICs')
plot_roe_state_1(state_j2_nodrift)

%% Overlay with numerical and STM
figure('Name','Numerical and STM overlay')
plot_roe_state_2(osc_qns_j2_roe_state,state_j2)
legend('Numerical','STM')

