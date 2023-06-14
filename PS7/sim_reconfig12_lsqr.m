%% Reconfiguration simulation from mode 1 to mode 2 using least squares
% AA 279D Problem Set 5
% Sydney Hsu and Pol Francesch
clear; clc; close all;

initial_conditions;

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
 
% Separate into desired reconfiguration ROE sets for in-plane and
% out-of-plane solutions
droe1to2_ip = droe1to2;
droe1to2_ip(5) = 0; droe1to2_ip(6) = 0;
droe1to2_oop = droe1to2;
droe1to2_oop(3) = 0; droe1to2_oop(4) = 0;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 Part b
% Least squares solution to reconfigure from C1 to C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Decoupled in-plane and out-of-plane 
% In-plane solution
N = 2;
desired_phase_ip = wrapTo2Pi(atan2(dey_2-dey_1,dex_2-dex_1)); % [from Chernick thesis]
%desired_phase_ip = wrapTo2Pi(-atan2(dex_2-dex_1,dey_2-dey_1)); % [from D'Amico]
k = 1; % natural number to select maneuver location
u_ks_ip = [desired_phase_ip, desired_phase_ip + k*pi]; % selected maneuver locations along an orbit
u_ks_ip_wrapped = wrapTo2Pi(u_ks_ip);
% change over full maneuver cycle
delta_us_ip = (u_ks_ip(end)-0)*ones(1,N); % location differences (mean arg. of lat) between maneuvers
t_man_ip = u_ks_ip / (2*pi) * T; % time of maneuvers

% Apply least-squares maneuver
[dV_lsqr_ip, dV_lsqr_err_ip] = least_squares_maneuvers(mu, J2, rE, N, u_ks_ip_wrapped, delta_us_ip, droe1to2_ip, oe_init_c);
dV_ip = reshape(dV_lsqr_ip,[3,N]);

% Out-of-plane solution
N = 1;
desired_phase_oop = wrapTo2Pi(atan2(diy_2-diy_1,dix_2-dix_1)); % [from Chernick thesis]
u_ks_oop = [desired_phase_oop]; % selected maneuver locations along an orbit
u_ks_oop_wrapped = wrapTo2Pi(u_ks_oop);
% change over full maneuver cycle
delta_us_oop = (u_ks_oop(end)-0)*ones(1,N); % location differences (mean arg. of lat) between maneuvers
%delta_us = zeros(1,N);
t_man_oop = u_ks_oop / (2*pi) * T; % time of maneuvers

% Apply least-squares maneuver
[dV_lsqr_oop, dV_lsqr_err_oop] = least_squares_maneuvers(mu, J2, rE, N, u_ks_oop_wrapped, delta_us_oop, droe1to2_oop, oe_init_c);
dV_oop = reshape(dV_lsqr_oop,[3,N]);

% Combine into maneuver scheduling queue
dV = [dV_oop, dV_ip];
t_man = [t_man_oop, t_man_ip];


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time span to validate over
orb_rev = 1.5;
% stepSize = T/100;
stepSize = 10.0;
% tspan = 0:60.0:3 * 86400.0;
tspan = 0:stepSize:orb_rev*T;
steps = length(tspan);

% Looping vars
roe = roe_init1;
oe_c = oe_init_c;
eci_c_t = rv_init_c; 
eci_d_t = rv_init_d;

% Storage variables
u   = zeros(steps,3);
eci_c = zeros(steps+1,6);
eci_d = zeros(steps+1,6);
qns_oe_c = zeros(steps+1,6);
qns_oe_d = zeros(steps+1,6);
qns_roe  = zeros(steps+1,6);

eci_c(1,:) = eci_c_t;
eci_d(1,:) = eci_d_t;
qns_oe_c(1,:) = qns_oe_init_c;
qns_oe_d(1,:) = qns_oe_init_d1;
qns_roe(1,:) = roe_init1;

% Looping
for j=1:steps
    % Check if we are close to a maneuver time
    if ~isempty(t_man) && abs(t_man(1) - tspan(j)) <= stepSize
        u(j,:) = dV(:,1)' / stepSize;

        t_man = t_man(2:end);
        dV    = dV(:,2:end);
    end
    
    % Apply computed control and timestep
    tspan_local = 0:stepSize/2:stepSize;
    [~,eci_c_tspan] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_local,eci_c_t,options);
    [~,eci_d_tspan] = ode45(@(t,z) dfq_ptb_sam(t,z,rE,mu,J2, u(j,:), eci_d_t),tspan_local,eci_d_t,options);

    eci_c_t = eci_c_tspan(end,:)'; % take the last step
    eci_d_t = eci_d_tspan(end,:)';

    [osc_sing_oe, mean_sing_oe, osc_oe,mean_oe,osc_roe,mean_roe] = rv2osc_and_mean(mu, rE, J2, 1, 0, eci_c_t, eci_d_t);

    % Save relevant variables
    eci_c(j+1,:) = eci_c_t;
    eci_d(j+1,:) = eci_d_t;
    qns_oe_c(j+1,:) = osc_oe(1,:);
    qns_oe_d(j+1,:) = osc_oe(2,:);
    qns_roe(j+1,:)  = osc_roe;
end

eci_c = eci_c(1:end-1,:);
eci_d = eci_d(1:end-1,:);
qns_oe_c = qns_oe_c(1:end-1,:); % osculating QNS
qns_oe_d = qns_oe_d(1:end-1,:);
qns_roe  = qns_roe(1:end-1,:);

% Convert position and velocity to RTN
[rho, rhodot] = absToRel(mu,tspan,eci_c,eci_d);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
earth_plot(rE,eci_c,eci_d); % Verify ECI propagator

% Plotting RTN
figure('Name', 'RTN')
plot_rtn(rho);

orbit_periods = tspan/T;
qns_labels = ["$a$ (m)", "$e_x$", "$e_y$", "$i$ (deg)",...
                "$\Omega$ (deg)", "$u$ (deg)"];
qns_roe_labels = ["$a\delta a$ (m)", "$a\delta \lambda$ (m)", "$a\delta e_x$ (m)", ...
               "$a\delta e_y$ (m)", "$a\delta i_x$ (m)", "$a\delta i_y$ (m)"];

% QNS absolute chief orbital elements
figure('Name','Chief QNS Absolute')
qns_oe_c_deg = qns_oe_c;
qns_oe_c_deg(:,4) = rad2deg(qns_oe_c(:,4));
qns_oe_c_deg(:,5) = rad2deg(qns_oe_c(:,5));
qns_oe_c_deg(:,6) = rad2deg(qns_oe_c(:,6));
plot_oe(orbit_periods,qns_oe_c_deg,qns_labels)
% QNS absolute deputy oes
qns_oe_d_deg = qns_oe_d;
qns_oe_d_deg(:,4) = rad2deg(qns_oe_d(:,4));
qns_oe_d_deg(:,5) = rad2deg(qns_oe_d(:,5));
qns_oe_d_deg(:,6) = rad2deg(qns_oe_d(:,6));
figure('Name','Deputy QNS Absolute')
plot_oe(orbit_periods,qns_oe_d_deg,qns_labels)

% QNS ROEs
qns_roe_norm = sma_c*qns_roe;
figure('Name','QNS ROE')
plot_oe(orbit_periods,qns_roe_norm,qns_roe_labels);
% figure('Name','ROE State Space')
% plot_roe_state(qns_roe_norm);

% Plot initial and final states in state space
figure('Name','ROE State Space Evolution');
roe_desired = roe_init2*sma_c;
plot_roe_state_w_init_and_final(qns_roe_norm,roe_init1*sma_c,roe_desired);

% Plot control acceleration inputs
figure('Name','Control Input')
subplot(1,2,1);
plot(orbit_periods,u(:,2));
xlabel('Orbit Periods')
ylabel('$u_T (m/s^2)$')
grid on;
subplot(1,2,2)
plot(orbit_periods,u(:,3));
xlabel('Orbit Periods')
ylabel('$u_N (m/s^2)$')
grid on;

% Plot cumulative delta-V over time
deltaVs = u*stepSize;
deltaVMag = vecnorm(deltaVs')';
deltaVcumulative = cumsum(deltaVMag);
figure('Name','Delta-V')
plot(orbit_periods,deltaVcumulative);
xlabel('Orbit Periods')
ylabel('Cumulative $\Delta v (m/s)$')
grid on;

% Plot control tracking error
roe_error = qns_roe_norm - roe_desired';
figure('Name','Control Tracking Error')
for j=1:length(roe_init2)
    plot(orbit_periods,roe_error(:,j)); hold on
end
xlabel('Orbit Periods')
ylabel('ROE error [m]')

qns_roe_labels = ["$a\delta a$", "$a\delta \lambda$", "$a\delta e_x$", ...
               "$a\delta e_y$", "$a\delta i_x$", "$a\delta i_y$"];
legend(qns_roe_labels,'Interpreter','latex');
grid on;


%% Reconfiguration accuracy
roe_desired = round(roe_desired,10);
% Achieved final ROE (at the end of propagation)
roe_achieved = qns_roe_norm(end,:);

% ROE final error in meters
roe_error_m = (roe_achieved'-roe_desired)

% Error between the achieved and desired (decimal)
rel_error = roe_error_m ./ roe_desired


%% Delta-V lower bounds
% In-plane maneuvers
delta_de = norm(droe1to2(3:4));
eta = sqrt(1-ecc_c^2);
deltaV_LB = n*sma_c*eta*delta_de/(2*eta^2)

dV_ip_total = sum(vecnorm(dV_ip,2,1))

% Out-of-plane lower bound
delta_di = norm(droe1to2(5:6));
deltaV_LB_oop = sma_c*n*(1-ecc_c)/eta*delta_di

dV_oop_total = sum(dV_oop)

