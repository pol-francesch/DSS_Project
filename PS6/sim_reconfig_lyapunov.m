%% Continuous control reconfiguration simulation from mode 1 to mode 2
% AA 279D Problem Set 5
% Sydney Hsu and Pol Francesch
clear; clc; close all;

initial_conditions;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desired change in ROEs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Setup reconfiguration: define initial orbit of chief, initial ROE, final ROE
% Final ROE conditions
% Mode 2: phase C2
dsma_2 = 0;
dlambda_2 = 0;
de_2 = 297 / sma_c;
di_2 = 260 / sma_c;
phase_diff_2 = 210; % designed phase difference % TN plane not relevant as it is safe in the other planes; along-track is the greatest uncertainty
theta_2 = theta_1; % phase angle of inclination

% Calculate deputy orbit parameters
phi_2 = theta_2 - deg2rad(phase_diff_2); % phase angle of eccentricity
dex_2 = de_2*cos(phi_2);  dey_2 = de_2*sin(phi_2);
dix_2 = di_2*cos(theta_2); diy_2 = di_2*sin(theta_2);
% Initial conditions - quasi-nonsingular ROEs
roe_init2 = [dsma_2;dlambda_2;dex_2;dey_2;dix_2;diy_2];
roe_init2_trunc = [dsma_2;dex_2;dey_2;dix_2;diy_2];
droe1to2 = roe_init2 - roe_init1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find optimal thrust using continous control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time span to validate over
orb_rev = 20;
stepSize = T/100;
tspan = 0:stepSize:T*orb_rev;
steps = length(tspan);

% Controller gain parameters
N = 20; % [Steindorf et al]
k_fdbk = 1000; % arbitrarily large scalar; TODO if too small, gain is too large and the orbit becomes unbounded (e>1)

% QNS ROE [dsma;dex;dey;dix;diy];
roe_init1_trunc = [roe_init1(1); roe_init1(3:6)];
roe = roe_init1_trunc;
oe_c = oe_init_c;
eci_c_t = rv_init_c; 
eci_d_t = rv_init_d;

% roe_init2_trunc = roe_init1_trunc;
% roe_init2       = roe_init1;

% Reference governor ROE to be updated each step
roe_refgov = roe_init2_trunc;

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

% Calculate delta-V lower bound
delta_de = norm(roe_init2_trunc(2:3)-roe_init1_trunc(2:3));
delta_di = norm(roe_init2_trunc(4:5)-roe_init1_trunc(4:5));
eta = sqrt(1-ecc_c^2);
deltaV_LB = (sma_c*n)/eta*(delta_de/2+(1-ecc_c)*delta_di);


for j=1:steps
    % Calculate plant & control matrix
    A = build_A(mu, oe_c);
    B = build_B(mu, oe_c);
    % [A, B] = get_reduced_model(oe_c);
    
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

%% Plot results
earth_plot(rE,eci_c,eci_d); % Verify ECI propagator

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
function A = build_A(mu, oe_c)
    sma_c = oe_c(1);
    n = sqrt(mu/sma_c^3); % mean motion

    % Keplerian
    A = zeros(6,6);
    A(2,1) = 1.5*n;
end

function B = build_B(mu, oe_c)
    % Unwrap state
    sma_c = oe_c(1); ecc_c = oe_c(2); inc_c = oe_c(3);
    raan_c = oe_c(4); w_c = oe_c(5); M_c = oe_c(6);

    ex_c = ecc_c * cos(w_c);
    ey_c = ecc_c * sin(w_c);

    n = sqrt(mu/sma_c^3); % mean motion

    ta_c = mean2true(M_c, ecc_c, 1e-12);

    % Get auxiliaries
    sfc = sin(ta_c);
    cfc = cos(ta_c);
    tic = tan(inc_c);
    swfc = sin(w_c + ta_c);
    cwfc = cos(w_c + ta_c);
    eta = sqrt(1 - ecc_c^2);

    % Get B
    B = zeros(5,2);

    B(1,1) = 2*ecc_c*sfc/eta;
    B(1,2) = 2*(1+ecc_c*cfc)/eta;
    B(2,1) = -2*eta^2/(1+ecc_c*cfc);
    B(3,1) = eta*swfc;
    B(3,2) = eta*((2+ecc_c*cfc)*cwfc+ex_c)/(1+ecc_c*cfc);
    B(3,3) = eta*ey_c*swfc/(tic*(1+ecc_c*cfc));
    B(4,1) = -eta*cwfc;
    B(4,2) = eta*((2+ecc_c*cfc)*swfc+ey_c)/(1+ecc_c*cfc);
    B(4,3) = -eta*ex_c*swfc/(tic*(1+ecc_c*cfc));
    B(5,3) = eta*cwfc/(1+ecc_c*cfc);
    B(6,3) = eta*swfc/(1+ecc_c*cfc);

    B = 1/(n*sma_c) * B;
end

function P = build_P(roe, roe_ref, oe_c, N, K)
    phi_ip = atan2(roe(3) - roe_ref(3),roe(2) - roe_ref(2)); % "Delta delta e": variation in e across the reconfig.
    phi_oop = atan2(roe(5) - roe_ref(5), roe(4) - roe_ref(4));
    
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