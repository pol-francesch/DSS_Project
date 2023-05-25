%% Reconfiguration simulation from mode 1 to mode 2 using D'Amico PhD
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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 Part b
% D'Amico PhD solution to reconfigure from C1 to C2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In-plane (single maneuver solution)
% de_1to2 = -norm(droe1to2(3:4)); % desired change in eccentricity
% if droe1to2(1) == 0
%     dv_t = 0;
%     dv_r = n*sma_c*de_1to2;
%     u_M_ip = - atan2(droe1to2(3),droe1to2(4));
% else
%     dv_t = n/2*sma_c*droe1to2(1);
%     dv_r = n*sma_c*sqrt(de_1to2^2+droe1to2(1)^2);
%     u_M_ip = atan2(dv_r,(2*dv_t)) - atan2(droe1to2(4),droe1to2(3));
% end
% 
% % Out-of-plane
% di_1to2 = sqrt(droe1to2(5)^2+droe1to2(6)^2);
% dv_n = n*sma_c*di_1to2;
% u_M_oop = atan2(droe1to2(6),droe1to2(5));
% 
% % t_man = [u_M_ip/(2*pi)*T,u_M_oop/(2*pi)*T+T]; % + T; % applied after 1 orbit
% % dV = [[dv_r;dv_t;0],[0;0; dv_n]];
% 
% t_man = [u_M_oop, u_M_ip]/(2*pi)*T;
% dV = [[0;0; dv_n],[dv_r;dv_t;0]];
% 
% dV_store = dV;
% t_man_store = t_man;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 Part b
% D'Amico PhD 2.4.6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Guidance - 2.4.4
delta_uM = pi; % spacing between consecutive burns; minimum delta-V cost solution by GVE
dinc_max = 0; % maximum allowed deviation in inclination
decc_max = 0; % maximum allowed deviation in eccentricity
dt = 3 * T; % maneuver cycle; time interval between the first pulses of 2 consecutive maneuver pairs

eta = sqrt(1 - ecc_c^2);
gamma = J2/2*(rE/sma_c)^2*1/eta^4;
phi_p = 3/2*gamma*(5*cos(inc_c)^2 - 1); % derivative of the rel. arg. of perigee w.r.t. mean arg. of latitude (D'Amico PhD eq. 2.30)
dphi = sign(phi_p)*asin(decc_max / norm(roe_init2(3:4))) + phi_p*delta_uM; % for 2-impulse control scheme

% Desired ecc/inc vectors after the in-plane maneuver execution
decc_man = [roe_init2(3)*cos(dphi) - roe_init2(4)*sin(dphi);
            roe_init2(3)*sin(dphi) + roe_init2(4)*cos(dphi)];
dinc_man = [roe_init2(5);
            roe_init2(6) - sign(roe_init1(5))*dinc_max];

% Desired sma after maneuvers
du_j2       = -12*gamma*sin(2*inc_c)*roe_init1(5)*n*dt;
dsma_man    = -pi/(2*n*dt - pi)*(3/2*norm(decc_man - roe_init1(3:4)) + roe_init1(1) - 4/(3*pi)*(roe_init1(2) - roe_init2(2) + du_j2));
% dsma_man = 0;

% Double maneuver solution with radial maneuvers - 2.4.6
% u_M1  = wrapTo2Pi(atan2(decc_man(1) - roe_init1(3), decc_man(2) - roe_init1(4)));
% u_M2  = delta_uM + u_M1;
% du_dv = -3*pi/4*(dsma_man-roe_init1(1) + norm(decc_man-roe_init1(3:4)));
% % du_T = 3*pi/8*norm(decc_man-roe_init1(3:4))+roe_init1(6);
% du_da = -3*pi/2*roe_init1(1);
% du_D = 0;
% du = roe_init1(2); % current rel. mean arg. of latitude right before the first maneuver
% dlambda_man = roe_init1(6) - du - du_da - du_j2 - du_D;
% %dlambda_man = 0;
% dv_r1 = n*sma_c/2 * ( -(dlambda_man - roe_init1(2))/2 + norm(decc_man - roe_init1(3:4)) );
% dv_r2 = n*sma_c/2 * ( -(dlambda_man - roe_init1(2))/2 - norm(decc_man - roe_init1(3:4)) );
% 
% % Compensate for semi-major axis changes
% dv_t1 = n*sma_c/4*(dsma_man-roe_init1(1));
% dv_t2 = n*sma_c/4*(dsma_man-roe_init1(1));
% 
% t_man = [u_M1, u_M2] / (2*pi)*T;
% dV = [ [dv_r1;dv_t1;0], [dv_r2;dv_t2;0] ];


% In-plane (double maneuver solution) - 2.4.4
u_M1  = wrapTo2Pi(atan2(decc_man(2) - roe_init1(4), decc_man(1) - roe_init1(3))); % desired change in phase of ecc.
u_M2  = delta_uM + u_M1;
dv_t1 = n*sma_c/4*((dsma_man - roe_init1(1)) + norm(decc_man - roe_init1(3:4)));
dv_t2 = n*sma_c/4*((dsma_man - roe_init1(1)) - norm(decc_man - roe_init1(3:4)));

% Out-of-plane
dv_n = n*sma_c*norm(dinc_man - roe_init1(5:6));
u_M_oop = atan2(dinc_man(2) - roe_init1(6), dinc_man(1) - roe_init1(5));

% Put maneuver schedule together
t_man = [u_M_oop, u_M1, u_M2]./(2*pi).*T;

dV = [[0;0; dv_n],...
      [0;dv_t1;0],...
      [0;dv_t2;0],...
      ];

dV_store = dV;
t_man_store = t_man;

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

    % Update looping vars
    roe = mean_roe';
    oe_c = mean_sing_oe(1,:);

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
deltaV_LB_ip = n*sma_c*eta*delta_de/(2*eta^2)

dV_ip_total = abs(dv_t1) + abs(dv_t2)

% Out-of-plane lower bound
delta_di = norm(droe1to2(5:6));
deltaV_LB_oop = sma_c*n*(1-ecc_c)/eta*delta_di

dV_oop_total = abs(dv_n)
