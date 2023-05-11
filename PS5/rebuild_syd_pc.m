%% AA 279D Problem Set 4
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
% Problem 1 - Part a
% Chief Initial Orbital Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial chief singular orbital parameters
sma_c  = 6892.927e3;      % semi-major axis [m]
ecc_c  = 1e-4;            % eccentricity
inc_c  = deg2rad(97.44);  % inclination [rad]
raan_c = deg2rad(270);    % RAAN [rad]
aop_c  = deg2rad(0);     % aop [rad]
ta_c   = deg2rad(0);    % true anomaly [rad]

oe_c = [sma_c, ecc_c, inc_c, raan_c, aop_c, ta_c];
oe_osc_j2_c = mean2osc(oe_c, 1);

% Chief QNS OE
qns_oe_c = singular2qns(oe_c);

% Initial position and velocity in inertial frame
[r_c,v_c] = singular_oe2rv(mu,oe_c);
[r_c_j2, v_c_j2] = singular_oe2rv(mu, oe_osc_j2_c);

% Orbital period
T = 2*pi*sqrt(sma_c^3/mu); % [sec]

% Verify singular oe & qns oe
oe_verify_rv  = rv2singular_oe(mu, r_c, v_c);
oe_verify_qns = qns2singular(qns_oe_c);
oe_verify_osc = osc2mean(oe_osc_j2_c);
oe_verify_osc_j2_rv = osc2mean(rv2singular_oe(mu, r_c_j2, v_c_j2),1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part b
% Deputy Initial Orbital Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Relative quasi non-singular oe [dimensions are m -> sma_c*dx]
dsma = 0;
dlambda = 100;
dex = 50;
dey = 100;
dix = 30;
diy = 200;

% Deputy QNS OE
qns_oe_d = qns_oe_c;
qns_oe_d(1) = qns_oe_d(1) + qns_oe_c(1)*dsma/qns_oe_c(1);
qns_oe_d(2) = qns_oe_d(2) + dex/qns_oe_c(1);
qns_oe_d(3) = qns_oe_d(3) + dey/qns_oe_c(1);
qns_oe_d(4) = qns_oe_d(4) + dix/qns_oe_c(1);
qns_oe_d(5) = qns_oe_d(5) + diy/qns_oe_c(1)/sin(qns_oe_c(4));
qns_oe_d(6) = qns_oe_d(6) + dlambda / qns_oe_c(1) - (qns_oe_d(5) - qns_oe_c(5))*cos(qns_oe_c(4));

% Deputy Singular OE
oe_d = qns2singular(qns_oe_d);
oe_osc_j2_d = mean2osc(oe_d,1);

% Initial position and velocity in inertial frame
[r_d,v_d] = singular_oe2rv(mu,oe_d);
[r_d_j2, v_d_j2] = singular_oe2rv(mu, oe_osc_j2_d);

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
z0_chief = [r_c;v_c];
z0_deputy = [r_d;v_d];

z0_chief_j2 = [r_c_j2; v_c_j2];
z0_deputy_j2 = [r_d_j2; v_d_j2];

% Numerical integration
[t,z_c_unptb] = ode45(@(t,z) dfq(t,z,rE,mu,0),tspan,z0_chief,options);
[~,z_d_unptb] = ode45(@(t,z) dfq(t,z,rE,mu,0),tspan,z0_deputy,options);

[~,z_c_j2] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan,z0_chief_j2,options);
[~,z_d_j2] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan,z0_deputy_j2,options);

% Verification
if false
    earth_plot(rE, z_c_unptb, z_d_unptb);
    earth_plot(rE, z_c_j2, z_d_j2);
end

%% To OE & ROE
steps = length(tspan);

% Third dimension is (chief, deputy)
osc_qns_unptb = zeros(steps, 6, 2);
mean_qns_unptb = zeros(steps, 6, 2);

osc_qns_unptb_roe = zeros(steps, 6);
mean_qns_unptb_roe = zeros(steps, 6);

osc_qns_j2 = zeros(steps, 6, 2);
mean_qns_j2 = zeros(steps, 6, 2);

osc_qns_j2_roe = zeros(steps, 6);
mean_qns_j2_roe = zeros(steps, 6);

for j=1:1:steps
    %%%%%%%%%%% Unperturbed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    chief = z_c_unptb(j,:); deputy = z_d_unptb(j, :);

    [osc_oe, mean_oe, osc_roe, mean_roe] = rv2osc_and_mean(mu, rE, J2, 0, 0, chief, deputy);
    
    osc_qns_unptb(j,:,:) = osc_oe';
    mean_qns_unptb(j,:,:) = mean_oe';
    osc_qns_unptb_roe(j,:) = osc_roe';
    mean_qns_unptb_roe(j,:) = mean_roe';

    %%%%%%%%%%% Perturbed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    chief = z_c_j2(j,:); deputy = z_d_j2(j, :);

    [osc_oe, mean_oe, osc_roe, mean_roe] = rv2osc_and_mean(mu, rE, J2, 1, 0, chief, deputy);
    
    osc_qns_j2(j,:,:) = osc_oe';
    mean_qns_j2(j,:,:) = mean_oe';
    osc_qns_j2_roe(j,:) = osc_roe';
    mean_qns_j2_roe(j,:) = mean_roe';
end

%% Plot OE & ROE
orbit_periods = t/T;
abs_ylabels = ["$a$ (m)", "$e_x$", "$e_y$", "$i$ (rad)",...
                "$\Omega$ (rad)", "$u$ (rad)"];
rel_ylabels = ["$\delta a$ (m)", "$\delta \lambda$", "$\delta e_x$", ...
               "$\delta e_y$", "$\delta i_x$", "$\delta i_y$"];

% Unperturbed absolute chief
plot_osc_mean(orbit_periods, mean_qns_unptb(:,:,1), osc_qns_unptb(:,:,1), abs_ylabels);

% Unperturbed absolute deputy
plot_osc_mean(orbit_periods, mean_qns_unptb(:,:,2), osc_qns_unptb(:,:,2), abs_ylabels);

% Unperturbed relative
plot_osc_mean(orbit_periods, mean_qns_unptb_roe, osc_qns_unptb_roe, rel_ylabels);

% Perturbed absolute chief
plot_osc_mean(orbit_periods, mean_qns_j2(:,:,1), osc_qns_j2(:,:,1), abs_ylabels);

% Perturbed absolute deputy
plot_osc_mean(orbit_periods, mean_qns_j2(:,:,2), osc_qns_j2(:,:,2), abs_ylabels);

% Perturbed relative
plot_osc_mean(orbit_periods, mean_qns_j2_roe, osc_qns_j2_roe, rel_ylabels);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part 4
% RTN frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get data
[rho2_unptb, ~] = absToRel(mu, t, z_c_unptb, z_d_unptb);
[rho2_j2, ~] = absToRel(mu, t, z_c_j2, z_d_j2);

% Plot
figure()
subplot(2,2,1); hold on;
plot(rho2_unptb(:,2),rho2_unptb(:,1), 'LineWidth',2); 
plot(rho2_j2(:,2),rho2_j2(:,1), '--'); 
title('TR'); xlabel('y (m)'); ylabel('x (m)'); axis equal; grid on
subplot(2,2,2); hold on;
plot(rho2_unptb(:,3),rho2_unptb(:,1), 'LineWidth',2); 
plot(rho2_j2(:,3),rho2_j2(:,1), '--'); 
title('NR'); xlabel('z (m)'); ylabel('x (m)'); axis equal; grid on
subplot(2,2,3); hold on;
plot(rho2_unptb(:,2),rho2_unptb(:,3), 'LineWidth',2); 
plot(rho2_j2(:,2),rho2_j2(:,3), '--'); 
title('TN'); xlabel('y (m)'); ylabel('z (m)'); axis equal; grid on
subplot(2,2,4); hold on; view(3);
plot3(rho2_unptb(:,1),rho2_unptb(:,2),rho2_unptb(:,3), 'LineWidth',2); 
plot3(rho2_j2(:,1), rho2_j2(:,2),rho2_j2(:,3), '--'); 
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); grid on
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

% Deputy QNS OE
qns_oe_d = qns_oe_c;
qns_oe_d(1) = qns_oe_d(1) + qns_oe_c(1)*dsma_new/qns_oe_c(1);
qns_oe_d(2) = qns_oe_d(2) + dex_new/qns_oe_c(1);
qns_oe_d(3) = qns_oe_d(3) + dey_new/qns_oe_c(1);
qns_oe_d(4) = qns_oe_d(4) + dix_new/qns_oe_c(1);
qns_oe_d(5) = qns_oe_d(5) + diy_new/qns_oe_c(1)/sin(qns_oe_c(4));
qns_oe_d(6) = qns_oe_d(6) + dlambda_new / qns_oe_c(1) - (qns_oe_d(5) - qns_oe_c(5))*cos(qns_oe_c(4));

% Deputy Singular OE
oe_d = qns2singular(qns_oe_d);
oe_osc_d = mean2osc(oe_d,1);

% Initial position and velocity in inertial frame
[r_d,v_d] = singular_oe2rv(mu,oe_osc_d);
z0_deputy = [r_d;v_d];

% Numerical integration
[t,z_c_j2_nodrift] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan,z0_chief_j2,options);
[t,z_d_j2_nodrift] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan,z0_deputy,options);

% Convert to osculating oe
osc_nodrift = zeros(steps,6,2);
mean_nodrift = zeros(steps,6,2);
osc_roe_nodrift = zeros(steps, 6);
mean_roe_nodrift = zeros(steps, 6);
for j=1:steps
    chief = z_c_j2_nodrift(j,:); deputy = z_d_j2_nodrift(j, :);

    [osc_oe, mean_oe, osc_roe, mean_roe] = rv2osc_and_mean(mu, rE, J2, 1, 0, chief, deputy);
    
    osc_nodrift(j,:,:) = osc_oe';
    mean_nodrift(j,:,:) = mean_oe';
    osc_roe_nodrift(j,:) = osc_roe';
    mean_roe_nodrift(j,:) = mean_roe';
end

% Visualization of relative drift
figure('Name','Drift Correction ROE')
orbit_periods = t/T;
for i = 1:1:6
    subplot(2,3,i); hold on;
    plot(orbit_periods, osc_roe_nodrift(:,i,1));
    xlabel('Orbital periods')
    ylabel(rel_ylabels(i));
    ytickformat('%.3f')
    grid on; hold off;
end

% State space plots
osc_state_roe = sma_c*osc_roe_nodrift;
plot_roe_state_1(osc_state_roe)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part 8
% J2 STM for ROE and comparison to numerical integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial conditions - quasi-nonsingular ROEs
state0 = [dsma;dlambda;dex;dey;dix;diy]; % units: [m]
state1_0 = [dsma_new;dlambda_new;dex_new;dey_new;dix_new;diy_new];

% Set up propagation
orb_rev = 15;
stepSize = T/1000;
tspan = 0:stepSize:T*orb_rev;
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances
state = zeros(length(tspan),6);
state1 = zeros(length(tspan),6);
for j=1:length(tspan)
    t = tspan(j);
    % Relative position and velocity
    state(j,:) = stm_qns_roe_j2(t,sma_c,ecc_c,inc_c,aop_c)*state0;
    state1(j,:) = stm_qns_roe_j2(t, sma_c, ecc_c, inc_c, aop_c)*state1_0;
end

figure('Name','J2 STM ROE State')
plot_roe_state_1(state)

figure('Name', 'J2 STM ROE State New ICs')
plot_roe_state_1(state1)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auxiliary Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function earth_plot(rE, chief, deputy)
    % Sanity check plot to verify orbit propagation
    [xE,yE,zE] = ellipsoid(0,0,0,rE,rE,rE,20);

    figure(); hold on; axis equal; grid on;
    set(gca,'DefaultLineLineWidth',1)
    plot3(chief(:,1),chief(:,2), chief(:,3), 'm');
    plot3(deputy(:,1),deputy(:,2), deputy(:,3), 'g');
    surface(xE,yE,zE, 'FaceColor','blue','FaceAlpha',.4,'EdgeColor','black','EdgeAlpha',0.5);
    title('Orbit around the Earth in ECI');
    xlabel('Position along I (m)');
    ylabel('Position along J (m)');
    zlabel('Position along K (m)')
    legend('Orbit 1', 'Orbit 2','Earth');
    view(3);
    hold off;
end

function plot_osc_mean(t, mean_elems, osc_elems, ylabels)
    % Plots osculating and mean elements
    % Works for either absolute or relative
    % 
    % Inputs:
    %   t           - time
    %   elems       - as many elements as you'd like
    % Outputs:
    
    figure();
    for i = 1:1:6
        subplot(2,3,i); hold on;
        plot(t, osc_elems(:,i));
        plot(t, mean_elems(:,i),'--');
        ylabel(ylabels(i));
        ytickformat('%.3f')
        grid on; hold off;
    end
    legend('Osculating','Mean')
end

% Plot overlaid unperturbed and J2 for dex-dey, dix-diy, dlambda-dsma
function plot_roe_state_2(unptb_roe,j2_roe)
    dims = size(unptb_roe);
    sets = dims(2) / 6;
    % Relative eccentricity
    subplot(3,1,1); hold on; grid on; axis equal;
    for i=1:sets
        plot(unptb_roe(:,3),unptb_roe(:,4)) %,'.','MarkerSize',10); hold on
        plot(j2_roe(:,3),j2_roe(:,4),'--') %,'.','MarkerSize',10);
        hold on
    end
    xlabel('$a \delta e_x $ (m)')
    ylabel('$a \delta e_y$ (m)')

    % Relative inclination
    subplot(3,1,2); hold on; grid on; axis equal;
    for i=1:sets
        plot(unptb_roe(:,5),unptb_roe(:,6),'.','MarkerSize',10); hold on
        plot(j2_roe(:,5),j2_roe(:,6),'.','MarkerSize',10);
        hold on
    end
    xlabel('$a \delta i_x $ (m)')
    ylabel('$a \delta i_y$ (m)')

    % Relative mean arg. of latitude and semi-major axis
    subplot(3,1,3); hold on; grid on; axis equal;
    for i=1:sets
        plot(unptb_roe(:,2),unptb_roe(:,1)); hold on
        plot(j2_roe(:,2),j2_roe(:,1),'--');
        hold on
    end
    xlabel('$a \delta \lambda $ (m)')
    ylabel('$a \delta a$ (m)')
end


% Plot dex-dey, dix-diy, dlambda-dsma
function plot_roe_state_1(roe_set)
    dims = size(roe_set);
    sets = dims(2) / 6;
    % Relative eccentricity
    subplot(3,1,1); hold on; grid on; axis equal;
    for i=1:sets
        u_i = roe_set(:,1+6*(i-1):6+6*(i-1));
        plot(u_i(:,3),u_i(:,4));
        hold on
    end
    xlabel('$a \delta e_x $ (m)')
    ylabel('$a \delta e_y$ (m)')

    % Relative inclination
    subplot(3,1,2); hold on; grid on; axis equal;
    for i=1:sets
        u_i = roe_set(:,1+6*(i-1):6+6*(i-1));
        plot(u_i(:,5),u_i(:,6),'.');
        hold on
    end
    xlabel('$a \delta i_x $ (m)')
    ylabel('$a \delta i_y$ (m)')

    % Relative mean arg. of latitude and semi-major axis
    subplot(3,1,3); hold on; grid on; axis equal;
    for i=1:sets
        u_i = roe_set(:,1+6*(i-1):6+6*(i-1));
        plot(u_i(:,2),u_i(:,1));
        hold on
    end
    xlabel('$a \delta \lambda $ (m)')
    ylabel('$a \delta a$ (m)')
end











