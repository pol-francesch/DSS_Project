%% AA 279D Problem Set 5
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

oe_init_c = [sma_c, ecc_c, inc_c, raan_c, aop_c, ta_c]; % combine the above into a single vector
qns_oe_init_c = singular2qns(oe_init_c)'; % chief QNS OE

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
roe_init1 = [dsma_1;dlambda_1;dex_1;dey_1;dix_1;diy_1]; % units: [m]

qns_oe_init_d1 = roe2qns(qns_oe_init_c,roe_init1); % QNS orbital elements of the deputy
oe_init_d1 = qns2singular(qns_oe_init_d1); % singular orbital elements of the deputy [sma, ecc, inc, raan, aop, ]
oe_init_d1_j2 = mean2osc(oe_init_d1, 1);
[r_init_d1,v_init_d1] = singular_oe2rv(mu,oe_init_d1); % deputy position and velocity in inertial frame, unptb
[r_init_d1_j2,v_init_d1_j2] = singular_oe2rv(mu,oe_init_d1_j2);

oe_init_c = qns2singular(qns_oe_init_c); % singular oe of chief
oe_init_c_j2 = mean2osc(oe_init_c, 1);
[r_init_c,v_init_c] = singular_oe2rv(mu,oe_init_c); % chief position and velocity in inertial frame, unptb
[r_init_c_j2,v_init_c_j2] = singular_oe2rv(mu,oe_init_c_j2);

% Orbital period of chief
T = 2*pi*sqrt(sma_c^3/mu); % [sec]

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 Orbit formation-keeping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% max. allowable relative position between satellites in RTN [m]
% ctrl_window_rtn = [250;250;250]; 

% Eccentricity control window
delta_psi_max = deg2rad(30);
e_vec = roe_init1(3:4);
e_max = norm(e_vec - e_vec*cos(delta_psi_max)); % unitless

% When we perform maneuvers, cannot exactly bounce back to the other side
% due to inexactness of dV
e_max_man = 0.9*e_max;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 Orbit formation-keeping
% Least squares maneuvers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orb_rev = 100;
days_elapsed = orb_rev*T/3600;
stepSize = T/1000;
tspan = 0:stepSize:T*orb_rev;

% Constant auxiliary values
eta = sqrt(1 - ecc_c^2);
gamma = J2/2 * (rE/sma_c)^2 * 1/eta^4;
phidot = 1.5*gamma*(5*cos(inc_c)^2 - 1);

ex_nom = roe_init1(3);
ey_nom = roe_init1(4);
e_norm = norm([ex_nom, ey_nom]);

n = sqrt(mu/sma_c^3);

% Set up propagation
steps = length(tspan);
% roe1 = zeros(steps,6);
man_times = [];
man_dVs   = [];

% If we perfom a maneuver, we will need to update the roe_0 that we use the
% STM off of. So we need to be able to keep a time since maneuver
% (local_tau) and a global time (tau).
tau = 0;
local_tau = 0;
roe_0 = roe_init1;

% for j=1:length(tspan)
while tau < orb_rev*T
    %tspan(j);
    % QNS ROE [dsma;dlambda;dex;dey;dix;diy];
    roe1 = stm_qns_roe_j2(local_tau,sma_c,ecc_c,inc_c,aop_c)*roe_0;

    % Check eccentricity control boundary
    if abs(norm(roe1(3:4) - [ex_nom; ey_nom])) < e_max
    % if abs(norm(roe1(3:4)) - e_norm) < e_max
        tau = tau + stepSize; 
        local_tau = local_tau + stepSize;
    else
        % Desired eccentricity after maneuver
        dphi = sign(phidot)*asin(e_max_man / e_norm);
        % dphi = sign(-phidot)*delta_psi_max*0.9;
        roe_ecc_man = [ex_nom*cos(dphi) - ey_nom*sin(dphi);
                       ex_nom*sin(dphi) + ey_nom*cos(dphi)];
        % roe_ecc_man = [ex_nom; ey_nom] / sin(delta_psi_max*0.9);
        
        % Calculate desired change in oe
        roe_man = roe1;
        % roe_man(3:4) = roe_init1(3:4);
        roe_man(3:4) = roe_ecc_man;
        droe = roe_man - roe1;

        % Calculate maneuver locations
        N = 2;

        if roe1(6) < pi
            % Currently on top of perifocal orbit. After periapsis before apo
            u_ks = [pi, 0];
            delta_uks = [pi - roe1(6), pi];
        else
            % After apo before peri
            u_ks = [0, pi];
            delta_uks = [2*pi - roe1(6), pi];
        end
        man_times_k = [tau + delta_uks(1)/n, tau + sum(delta_uks)/n];
        
        % Calculate maneuvers
        [dV, error] = least_squares_maneuvers(mu, J2, rE, N,u_ks,delta_uks, droe, oe_init_c);
        droe_actual = error + droe;

        % Save maneuvers
        man_times = [man_times, man_times_k];
        man_dVs   = [man_dVs, dV(1:3), dV(4:6)];

        % Update ROE and time
        tau = man_times_k(2);
        local_tau = 0;

        roe1 = droe_actual + roe1;
        roe_0 = roe1;
    end
end

sum(abs(man_dVs))

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 Orbit formation-keeping
% Validation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% orb_revs = 10;
% % days_elapsed = orb_revs*T/3600;
% stepSize = T/1000;
% tspan = 0:stepSize:T*orb_revs;
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances
rv_init_c = [r_init_c_j2,v_init_c_j2];
rv_init_d = [r_init_d1_j2,v_init_d1_j2];

% Propagate the chief and deputy orbits
[t, state_eci_c] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan,rv_init_c,options);
% [t,state_eci_c] = ode45(@(t,z) dfq_ptb(t,z,rE,mu,J2,0,0,0),tspan,rv_init_c,options);
% [~,state_eci_d] = ode45(@(t,z) dfq_ptb(t,z,rE,mu,J2,stepSize,man_times,man_dVs),tspan,rv_init_d,options);
% [~,state_eci_d] = ode45(@(t,z) dfq_ptb(t,z,rE,mu,J2,0,0,0),tspan,rv_init_d,options);

% % Initialize state arrays
state_eci_d = []; % zeros(steps,6);

% Manually propagate to each maneuver time
t0 = 0;
r0_d = rv_init_d;
true_tspan = [];

for j=1:length(man_times)
    % Propagate to maneuver
    tspan_to_man = t0:stepSize:man_times(j);
    [~, deputy] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_to_man,r0_d,options);

    % Perform maneuver
    post_man = deputy(end,:) + [0,0,0,man_dVs(:,j)'];

    % Add to common state
    if tspan_to_man(end) == man_times(j)
        deputy = deputy(1:end-1,:);
        tspan_to_man = tspan_to_man(1:end-1);
    end
    state_eci_d = [state_eci_d; deputy; post_man];
    true_tspan  = [true_tspan, tspan_to_man, man_times(j)];

    r0_d = post_man;
    t0   = man_times(j) + stepSize;
end

% tspan = true_tspan;
state_eci_d = interp1(true_tspan, state_eci_d, tspan);

earth_plot(rE,state_eci_c,state_eci_d);

%% Osculating orbital elements
osc_qns_c = zeros(steps,6);
osc_sing_c = zeros(steps,6);
osc_qns_d = zeros(steps,6);
osc_qns_roe = zeros(steps,6);

for j=1:steps
    oe_c_sing = rv2singular_oe(mu,state_eci_c(j,1:3),state_eci_c(j,4:6)); % pos/vel to singular oe
    osc_sing_c(j,:) = oe_c_sing;
    qns_c = singular2qns(oe_c_sing); % singular oe to QNS oe;
    osc_qns_c(j,:) = qns_c;
    oe_d_sing = rv2singular_oe(mu,state_eci_d(j,1:3),state_eci_d(j,4:6)); % pos/vel to singular oe
    qns_d = singular2qns(oe_d_sing); % singular oe to QNS oe
    osc_qns_d(j,:) = qns_d;
    osc_qns_roe(j,:) = oes2roe(qns_c,qns_d);
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
plot_roe_state_w_init_and_control_bounds(osc_qns_roe_norm, roe_init1*sma_c, e_max*sma_c);
% figure('Name','Singular OE Chief')
% plot_oe(orbit_periods,osc_sing_c,qns_labels);

% figure('Name','Initial State Space')
% osc_qns_roe_norm = sma_c*osc_qns_roe;
% plot_roe_state(osc_qns_roe_norm(1:50,:));
% 
% figure('Name','Final State Space')
% osc_qns_roe_norm = sma_c*osc_qns_roe;
% plot_roe_state(osc_qns_roe_norm(end-50:end,:));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 
% Helper fcns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_roe_state_w_init_and_control_bounds(roe_set, roe_init, ecc_max)
    % Eccentricity controller bounds
    th = 0:pi/50:2*pi;
    x_circ = ecc_max*cos(th) + roe_init(3);
    y_circ = ecc_max*sin(th) + roe_init(4);

    % Relative eccentricity
    subplot(3,1,1); hold on; grid on; axis equal;
    plot(roe_set(:,3),roe_set(:,4),'.','MarkerSize',10); hold on
    plot(roe_init(3), roe_init(4), '.', 'MarkerSize', 20);
    plot(x_circ, y_circ, 'r--');
    hold on
    xlabel('$a \delta e_x $ (m)')
    ylabel('$a \delta e_y$ (m)')

    % Relative inclination
    subplot(3,1,2); hold on; grid on; axis equal;
    plot(roe_set(:,5),roe_set(:,6),'.','MarkerSize',5); hold on
    plot(roe_init(5), roe_init(6), '.', 'MarkerSize', 10);
    hold on
    xlabel('$a \delta i_x $ (m)')
    ylabel('$a \delta i_y$ (m)')

    % Relative mean arg. of latitude and semi-major axis
    subplot(3,1,3); hold on; grid on; axis equal;
    plot(roe_set(:,2),roe_set(:,1)); hold on
    plot(roe_init(2), roe_init(1), '.', 'MarkerSize', 20);
    hold on
    xlabel('$a \delta \lambda $ (m)')
    ylabel('$a \delta a$ (m)')
    legend('Time history', 'Initial')
end
