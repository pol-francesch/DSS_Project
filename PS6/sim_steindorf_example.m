clear; clc; close all;

%% Preliminary setup
% Add function paths
addpath('.\plots');
addpath('.\mean_osc');
addpath('.\mappings');
addpath('.\dfqs');

% Constants
J2 = 0.0010826358191967; % Earth's J2 coefficient
mu = 3.986004415e14; % Earth gravitational parameter [m^3/s^2]
rE = 6.378136300e6; % Earth radius [m]

% Plotting
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesFontSize',18);

% Integration
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions - Chief
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial chief absolute singular orbital parameters
sma_c  = 6892e3;      % semi-major axis [m] %sma_c = 10000e3;
ecc_c  = 0.0001384;            % eccentricity component in I
inc_c  = deg2rad(97.44);  % inclination [rad]
raan_c = deg2rad(266.1539);    % RAAN [rad]
aop_c  = deg2rad(89.1198);      % aop [rad]
M_c    = deg2rad(45.88);      % mean anomaly [rad]

oe_init_c = [sma_c, ecc_c, inc_c, raan_c, aop_c, M_c]; % combine the above into a single vector
qns_oe_init_c = singular2qns(oe_init_c)'; % chief QNS OE (non-dimensionalized)
oe_init_c_j2 = mean2osc(oe_init_c, 1);

% Position and velocity
[r_init_c,v_init_c] = singular_oe2rv(mu,oe_init_c_j2); % chief position and velocity in inertial frame, J2-perturbed
rv_init_c = [r_init_c;v_init_c];

% Helpful parameters
T = 2*pi*sqrt(sma_c^3/mu); % [sec]
n = sqrt(mu/sma_c^3); % mean motion

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions - ROEs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
roe_init1_norm = [0;712;400;400;0;350];
roe_init2_norm = [0;212;200;200;0;200];
% roe_init2_norm = [0; 0; 0; 0; 0; 0];

delta_roe_norm = roe_init2_norm - roe_init1_norm;

roe_init1 = roe_init1_norm/sma_c;
roe_init2 = roe_init2_norm/sma_c;

roe_init2_trunc = [roe_init2(1); roe_init2(3:6)];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions - Deputy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial pos/vel of deputy
oe_init_d1 = roe2singular(oe_init_c,roe_init1);
qns_oe_init_d1 = singular2qns(oe_init_d1);

oe_init_d1_j2 = mean2osc(oe_init_d1, 1);

% Position and velocity
[r_init_d,v_init_d] = singular_oe2rv(mu,oe_init_d1_j2); % deputy position and velocity in inertial frame, J2 perturbed
rv_init_d = [r_init_d;v_init_d];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find optimal thrust using continous control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time span to validate over
orb_rev = 30;
stepSize = T/100;
tspan = 0:stepSize:T*orb_rev;
steps = length(tspan);

% Hyperparameters
N = 20; % [Steindorf et al]
k_fdbk = 1000; % arbitrarily large scalar; TODO if too small, gain is too large and the orbit becomes unbounded (e>1)

% QNS ROE [dsma;dex;dey;dix;diy];
roe_init1_trunc = [roe_init1(1); roe_init1(3:6)];
roe = roe_init1_trunc;
oe_c = oe_init_c;
eci_c_t = rv_init_c; 
eci_d_t = rv_init_d;

% Reference governor ROE to be updated each step
roe_refgov = roe_init2_trunc;

% Storage variables
u   = zeros(steps,2);
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

for j=1:steps
    % Calculate plant & control matrix
    [A, B] = get_reduced_model(oe_c);
    
    % Optimal control locations (update at each time step)
    P = build_P(roe, roe_refgov, oe_c, N, k_fdbk);

    % Compute control inputs
    tracking_error = roe - roe_refgov; % "Delta delta alpha"
    u(j,:) = -pinv(B)*(A*tracking_error + P*tracking_error);

    % Update chief & deputy
    tspan_local = 0:stepSize/2:stepSize;
    [~,eci_c_t] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_local,eci_c_t,options);
    % [~,eci_d_t] = ode45(@(t,z) dfq_ptb(t,z,rE,mu,J2, u(j,:), oe_c),tspan_local,eci_d_t,options);
    [~,eci_d_t] = ode45(@(t,z) dfq_ptb_works(t,z,rE,mu,J2, u(j,:), eci_d_t),tspan_local,eci_d_t,options);
    % [~,eci_d_t] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_local,eci_d_t,options);
    % [~,eci_d_t] = ode45(@(t,z) dfq_ptb(t,z,rE,mu,J2, [0,0], oe_c),tspan_local,eci_d_t,options);

    eci_c_t = eci_c_t(end,:)'; % take the last step
    eci_d_t = eci_d_t(end,:)';
    
    [osc_sing_oe, mean_sing_oe, osc_oe,mean_oe,osc_roe,mean_roe] = rv2osc_and_mean(mu, rE, J2, 1, 0, eci_c_t, eci_d_t);
    
    % Update looping vars
    roe = [mean_roe(1), mean_roe(3:6)]';
    oe_c = mean_oe(1,:);

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

%% Plot results
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
plot_roe_state_w_init_and_final(qns_roe_norm,roe_init1*sma_c,roe_init2*sma_c);

% Plot control acceleration inputs
figure('Name','Control Input')
subplot(1,2,1);
plot(orbit_periods,u(:,1));
xlabel('Orbit Periods')
ylabel('$u_T (m/s^2)$')
grid on;
subplot(1,2,2)
plot(orbit_periods,u(:,2));
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

% Calculate delta-V lower bound
delta_de = norm(roe_init2_trunc(2:3)-roe_init1_trunc(2:3));
delta_di = norm(roe_init2_trunc(4:5)-roe_init1_trunc(4:5));
eta = sqrt(1-ecc_c^2);
deltaV_LB = (sma_c*n)/eta*(delta_de/2+(1-ecc_c)*delta_di);

% Plot control tracking error
roe_error = qns_roe_norm - sma_c*roe_init2';
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


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = build_P(roe, roe_ref, oe_c, N, K)
    droe = roe - roe_ref;
    phi_ip = atan2(droe(3),droe(2)); % "Delta delta e": variation in e across the reconfig.
    phi_oop = atan2(droe(5), droe(4));
    
    % Current chief orbit
    u_M = oe_c(5) + oe_c(6); % mean arg. of latitude; aop + M

    % Enforce optimality
    J = u_M - phi_ip;
    H = u_M - phi_oop;
    P = 1/K*[cos(J)^N, 0, 0, 0, 0;
                  0, cos(J)^N, 0, 0, 0;
                  0, 0, cos(J)^N, 0, 0;
                  0, 0, 0, cos(H)^N, 0;
                  0, 0, 0, 0, cos(H)^N];
end

