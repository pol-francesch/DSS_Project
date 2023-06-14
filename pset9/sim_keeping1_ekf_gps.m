%% Formation keeping simulation for mode 1 using D'Amico's thesis
% AA 279D Problem Set 5
% Sydney Hsu and Pol Francesch
clear; clc; close all;

initial_conditions_pset9;
initial_conditions_pset9_gps;

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

% e_max = abs(norm(e_vec)*sin(delta_psi_max / sign(phi_p)));

% When we perform maneuvers, cannot exactly bounce back to the other side
% due to inexactness of dV
e_max_man = 0.9*e_max;

% Guidance - 2.4.4
delta_uM = pi;
dinc_max = 0;
decc_max = e_max_man;
dt = 50 * T;

eta = sqrt(1 - ecc_c^2);
gamma = J2/2*(rE/sma_c)^2*1/eta^4;
phi_p = 3/2*gamma*(5*cos(inc_c)^2 - 1);
% dphi = sign(-phi_p)*asin(decc_max / norm(roe_init1(3:4))) + phi_p*delta_uM;
dphi = asin(decc_max / norm(roe_init1(3:4))) + phi_p*delta_uM;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize EKF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(1); % initialize random seed

% Variance for ROE (m)
sigma = 10 * eye(6);

% Vary initial values by variance
qns_roe_ekf_t = (sqrtm(sigma)*randn(6,1) +  (roe_init1*sma_c))/sma_c;

% Matrices
P  = sigma; % estimation covariance
Q  = sigma / 1000; % process covariance
R  = sigma; % measurement covariance

% Measurement model
H = eye(6);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time span to validate over
orb_rev = 150;
stepSize = 10.0;
tspan = 0:stepSize:orb_rev*T;
steps = length(tspan);
uspan = zeros(steps,1); % mean arg. of latitude history (rad)

% Looping vars
roe = roe_init1; % quasi-nonsingular ROEs
eci_c_t = rv_init_c; 
eci_d_t = rv_init_d;

% Storage variables
u   = zeros(steps,3);
eci_c = zeros(steps+1,6);
eci_d = zeros(steps+1,6);
qns_oe_c = zeros(steps+1,6);
qns_oe_d = zeros(steps+1,6);
qns_osc_roe  = zeros(steps+1,6);

eci_c(1,:) = eci_c_t;
eci_d(1,:) = eci_d_t;
qns_oe_c(1,:) = qns_oe_init_c;
qns_oe_d(1,:) = qns_oe_init_d1;
qns_osc_roe(1,:) = roe_init1;

% Data storage for EKF
qns_roe_ekf  = zeros(steps+1,6);
qns_roe_ekf(1,:) = qns_roe_ekf_t;
% Store pre-fit and post-fit residuals
prefit_resids = zeros(steps+1,6);
postfit_resids = zeros(steps+1,6);
% Store the state error (covariance P)
P_diags = zeros(steps+1,6);

% u_tol = meanArgLatAtTime(stepSize,0,0,n,aop_c); % mean arg. of latitude tolerance

% Impulsive control
dV_queue = [];
t_man = [];
% u_man_queue = [];
u_true = 0;
bound_violations = 0;

% STM - invariant
sma0 = oe_init_c(1); ecc0 = oe_init_c(2); inc0 = oe_init_c(3); aop0 = oe_init_c(5);
stm = stm_qns_roe_j2(stepSize,sma0,ecc0,inc0,aop0);
B = impulsive_ctrl_input(n,eta,oe_init_c);

% Looping
for j=1:steps
    t = tspan(j);
    % Compute phase angle at current time step
    u_t = meanArgLatAtTime(t,0,0,n,aop_c); % mean arg. of latitude at current time
    
    % Check schedule to see if we are close to a maneuver time / location
%     if ~isempty(u_man_queue) && abs(u_man_queue(1) - u_t) <= u_tol
    if ~isempty(t_man) && abs(t_man(1) - tspan(j)) <= stepSize
        % IP maneuvers seem to be sensitive to timing inaccuracies
        if size(dV_queue,2) == 3
            dv_t1 = n*sma_c/4*((dsma_man - roe(1)) + norm(decc_man - roe(3:4)));
            dv_t2 = n*sma_c/4*((dsma_man - roe(1)) - norm(decc_man - roe(3:4)));
            dV_queue(:,1) = [0; dv_t1; 0];
            dV_queue(:,2) = [0; dv_t2; 0];
        end
        % OOP maneuvers aren't that sensitive to timing inaccuracies
        u(j,:) = dV_queue(:,1)' / stepSize;

        % Remove planned maneuver from the front of the queue
%         u_man_queue = u_man_queue(2:end); 
        t_man = t_man(2:end);
        dV_queue    = dV_queue(:,2:end);

    end

    % Check if we are out of control window
    % Only consider this if we have no maneuvers in buffer
    % We only perform maneuvers within one/two orbit so this is good enough
%     if isempty(u_man_queue) && abs(norm(roe(3:4) - roe_init1(3:4))) > e_max
    if abs(norm(roe(3:4) - roe_init1(3:4))) > e_max && isempty(t_man)
%     if abs(norm(roe(3:4) - roe_init1(3:4))) > e_max_man && isempty(t_man)
        % Guidance - 2.4.4
        decc_man = [roe_init1(3)*cos(dphi) - roe_init1(4)*sin(dphi);
                    roe_init1(3)*sin(dphi) + roe_init1(4)*cos(dphi)];
        dinc_man = [roe_init1(5);
                    roe_init1(6) - sign(roe(5))*dinc_max];
                
        du_j2       = -12*gamma*sin(2*inc_c)*roe(5)*n*dt;
        dsma_man    = -pi/(2*n*dt - pi)*(3*e_max_man + roe(1) - 4/(3*pi)*(roe(2) - roe_init1(2) + du_j2));

        % In-plane (double maneuver solution) - 2.4.4
        dv_t1 = n*sma_c/4*((dsma_man - roe(1)) + norm(decc_man - roe(3:4)));
        dv_t2 = n*sma_c/4*((dsma_man - roe(1)) - norm(decc_man - roe(3:4)));
        u_M1  = wrapTo2Pi(atan2(decc_man(2) - roe(4), decc_man(1) - roe(3)));
        u_M2  = delta_uM + u_M1;
        
        % Out-of-plane
        dv_n = n*sma_c*norm(dinc_man - roe(5:6));
        u_M_oop = wrapTo2Pi(atan2(dinc_man(2) - roe(6), dinc_man(1) - roe(5)));

        % Put it together
        dV_man_curr = [[0;dv_t1;0],...
                       [0;dv_t2;0],...
                       [0;0; dv_n]];
        % Normal maneuver times + current time + time to get to periapsis
        t_man_curr  = [u_M1/(2*pi)*T,...
                       u_M2/(2*pi)*T,...
                       u_M_oop/(2*pi)*T+T] + tspan(j) + (T - qns_oe_d(j,6) / (2*pi)*T);
        % Normal maneuver loc + current loc + angle remaining to get to periapsis
        % 2-orbit offset added so that oop maneuver occurs after both in-plane
%         u_man_curr  = [u_M1,u_M2,u_M_oop+4*pi] + u_t + (2*pi - qns_oe_d(j,6));
        
        % Ensure OOP comes last
        while t_man_curr(3) < t_man_curr(2)
            t_man_curr(3) = t_man_curr(3) + T;
        end

        % Store in buffer
        dV_queue = dV_man_curr;
        t_man = t_man_curr;
%         u_man_queue = u_man_curr;

        % Count how many maneuvers were needed
        bound_violations = bound_violations + 1;

    end
    
    % Apply computed control and timestep
    tspan_local = 0:stepSize/2:stepSize;
    [~,eci_c_tspan] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_local,eci_c_t,options);
    [~,eci_d_tspan] = ode45(@(t,z) dfq_ptb(t,z,rE,mu,J2, u(j,:), eci_d_t),tspan_local,eci_d_t,options);

    eci_c_t = eci_c_tspan(end,:)'; % take the last step
    eci_d_t = eci_d_tspan(end,:)';

    [osc_sing_oe, mean_sing_oe, osc_oe,mean_oe,osc_roe,mean_roe] = rv2osc_and_mean(mu, rE, J2, 1, 0, eci_c_t, eci_d_t);

    %%%% MEASUREMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    roe_meas = (sqrtm(sigma)*randn(6,1) +  (mean_roe'*sma_c))/sma_c; % corrupt the ROEs 

    %%%% EKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time update
    % Unpack mean singular orbital elements to determine the ROE state
    qns_roe_minus = stm * qns_roe_ekf_t + B*u(j,:)'*stepSize; % control inputs accounted for in prediction B*u
    P_minus = stm * P * stm' + Q;

    % Measurement update
    K = P_minus * H' / (H * P_minus * H' + R); % 
    prefit_residual = roe_meas - H*qns_roe_minus;
    qns_roe_ekf_t = qns_roe_minus + K*(prefit_residual); % QNS state at time t
    P = (eye(6) - K*H) * P_minus * (eye(6) - K*H)' + K*R*K'; % estimate covariance
    
    % Update looping vars
%     roe = mean_roe'; % zero noise
%     roe = roe_meas; % noisy measurements
    roe = qns_roe_ekf_t; % filtered measurements

    % Save relevant variables
    qns_roe_ekf(j+1,:) = qns_roe_ekf_t;
    P_diags(j+1,:) = diag(P);
    prefit_resids(j+1,:) = prefit_residual;
    postfit_resids(j+1,:) = roe_meas - H*qns_roe_ekf_t;
    uspan(j) = u_t; % store for reference
    eci_c(j+1,:) = eci_c_t;
    eci_d(j+1,:) = eci_d_t;
    qns_oe_c(j+1,:) = osc_oe(:,1);
    qns_oe_d(j+1,:) = osc_oe(:,2);
    qns_osc_roe(j+1,:)  = osc_roe; % ground truth

    u_true = u_true + (qns_oe_d(j+1,6) - qns_oe_d(j,6));
end

% Post-processing
eci_c = eci_c(1:end-1,:);
eci_d = eci_d(1:end-1,:);
qns_oe_c = qns_oe_c(1:end-1,:); % osculating QNS
qns_oe_d = qns_oe_d(1:end-1,:);
qns_osc_roe  = qns_osc_roe(1:end-1,:);

% Remove the very last row of data
qns_roe_ekf = qns_roe_ekf(1:end-1,:);
% Remove the first row of data
P_diags = P_diags(2:end,:);
prefit_resids = prefit_resids(2:end,:);
postfit_resids = postfit_resids(2:end,:);

% Convert position and velocity to RTN
[rho, rhodot] = absToRel(mu,tspan,eci_c,eci_d);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
earth_plot(rE,eci_c,eci_d); % Verify ECI propagator

% Plotting RTN
figure('Name', 'RTN')
plot_rtn(rho,'-',1);

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

%% QNS ROEs
qns_osc_roe_norm = sma_c*qns_osc_roe;
figure('Name','Osculating QNS ROE')
plot_oe(orbit_periods,qns_osc_roe_norm,qns_roe_labels);
% figure('Name','ROE State Space')
% plot_roe_state(qns_roe_norm);

% EKF
qns_roe_ekf_norm = qns_roe_ekf*sma_c;
figure('Name','QNS ROE EKF')
plot_oe(orbit_periods,qns_roe_ekf_norm,qns_roe_labels);

% Plot initial and final states in state space
roe_desired = roe_init1*sma_c;
% figure('Name','ROE State Space Evolution');
% plot_roe_state_w_init_and_final(qns_osc_roe_norm,roe_init1*sma_c,roe_desired);

%% Plot control acceleration inputs
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

%% Plot control tracking error
ctrl_error = qns_osc_roe_norm - roe_desired';
figure('Name','Control Tracking Error')
for j=1:length(roe_init1)
    plot(orbit_periods,ctrl_error(:,j)); hold on
end
xlabel('Orbit Periods')
ylabel('ROE error [m]')

qns_roe_labels = ["$a\delta a$", "$a\delta \lambda$", "$a\delta e_x$", ...
               "$a\delta e_y$", "$a\delta i_x$", "$a\delta i_y$"];
legend(qns_roe_labels,'Interpreter','latex');
grid on;

% Plot state space with control window
figure('Name', 'ROE State Space')
plot_roe_state_w_init_and_control_bounds(qns_osc_roe_norm, roe_init1*sma_c, e_max*sma_c);

%% Control tracking error steady-state statistics
% Steady-state defined as the last 10 orbits
ss_idx = find(orbit_periods>140,1);
ss_ctrl_data = ctrl_error(ss_idx:end,:);
% True mean/error (time series expectation)
ctrl_err_true_mean = mean(ss_ctrl_data);
ctrl_err_true_std = std(ss_ctrl_data);
% "formal"
% ctrl_error_mean_f = ;

% Map the control tracking error to RTN
[r_des_ECI,v_des_ECI] = singular_oe2rv(mu,roe2singular(oe_init_c,roe_init1));
[rtn_des_pos,rtn_des_vel] = absToRelState(mu,[r_init_c;v_init_d],[r_des_ECI;v_des_ECI]);
rtn_desired = [rtn_des_pos,rtn_des_vel];
ctrl_data_rtn = zeros(steps,6);
for j=1:steps
    % Convert ROEs to singular, singular to ECI pos/vel
    [r_ECI,v_ECI] = singular_oe2rv(mu,roe2singular(oe_init_c,qns_osc_roe(j,:)));
    % Convert ECI to RTN pos/vel
    [r_RTN,v_RTN] = absToRelState(mu,[r_init_c;v_init_c],[r_ECI;v_ECI]);
    ctrl_data_rtn(j,:) = [r_RTN,v_RTN];
end
ctrl_err_rtn = rtn_desired - ctrl_data_rtn; % relative position error

figure('Name','RTN Control Error'); hold on
rtn_labels = ["R (m)","T (m)","N (m)"];
plot_rtn_time_series(orbit_periods,ctrl_err_rtn(:,1:3),rtn_labels);

% Statistics in RTN
ctrl_err_rtn_mean = mean(ctrl_err_rtn);
ctrl_err_rtn_std = std(ctrl_err_rtn);

% Control accuracy in steady-state
ss_ctrl_error_rtn = ctrl_err_rtn(ss_idx:end,:);
ctrl_rms_3d = rms(vecnorm(ss_ctrl_error_rtn(:,1:3),2,2)) % root-mean squared 3D positional navigation accuracy [meters]


%% NAVIGATION PLOTS %%
% Plot navigation/estimation error with standard deviations/covariance bounds, split by ROE
nav_error = (qns_osc_roe - qns_roe_ekf) * sma_c;
P_3std_dev = 3*sqrt(P_diags);
figure('Name', 'Estimation error with bounds'); hold on;
plot_oe_error(orbit_periods,nav_error,P_3std_dev,qns_roe_labels)
legend('\pm3\sigma','State error')

% True navigation error statistics (computed from data)
% Time series expectation
ss_nav_data = nav_error(ss_idx:end,:);
nav_err_true_mean = mean(ss_nav_data);
nav_err_true_std = std(ss_nav_data); % conservative 

% TODO: formal navigation error statistics
nav_err_mean_formal = zeros(1,6); % mean error = 0 because of zero-mean Gaussian noise [NASA best practices]
nav_err_std_formal = sqrt(P_diags(end,:));

% Map the navigation error to RTN
[r_des_ECI,v_des_ECI] = singular_oe2rv(mu,roe2singular(oe_init_c,roe_init1));
rtn_desired = [r_des_ECI;v_des_ECI];
nav_data_rtn = zeros(steps,6);
for j=1:steps
    % Convert ROEs to singular, singular to ECI pos/vel
    [r_ECI,v_ECI] = singular_oe2rv(mu,roe2singular(oe_init_c,qns_roe_ekf(j,:))); 
    % Convert to RTN pos/vel
    [r_RTN,v_RTN] = absToRelState(mu,[r_init_c;v_init_c],[r_ECI;v_ECI]);
    nav_data_rtn(j,:) = [r_RTN,v_RTN];
end
nav_err_rtn = ctrl_data_rtn - nav_data_rtn; % relative navigation error

figure('Name','RTN Navigation Error'); hold on
plot_rtn_time_series(orbit_periods,nav_err_rtn(:,1:3),rtn_labels);

% Navigation accuracy in steady-state
ss_nav_error_rtn = nav_err_rtn(ss_idx:end,:);
nav_err_rtn_mean = mean(ss_nav_error_rtn);
nav_err_rtn_std = std(ss_nav_error_rtn);
nav_rms_3d = rms(vecnorm(ss_nav_error_rtn(:,1:3),2,2)) % root-mean squared 3D positional navigation accuracy [meters]

%% Zoom into the steady-state region
% figure('Name','Steady-State')
% idx = 1000; % 172
% % plot_oe(orbit_periods(1:idx),P_std_dev(1:idx,:),qns_roe_labels)
% plot_oe_error(orbit_periods(1:idx),nav_error(1:idx,:),P_std_dev(1:idx,:),qns_roe_labels)
% legend('\pm 3\sigma','State error')


%% Residuals
figure('Name','Residuals')
dims1 = size(prefit_resids);

% Normalize to sma
prefit_resids_norm = prefit_resids*sma_c;
postfit_resids_norm = postfit_resids*sma_c;

for i = 1:1:dims1(2)
    subplot(2,3,i); hold on;
    scatter(orbit_periods(ss_idx:end), prefit_resids_norm((ss_idx:end),i),3,'filled');
    scatter(orbit_periods(ss_idx:end), postfit_resids_norm((ss_idx:end),i),3,'filled');
    %xticks(1:1:steps)
%     xticks(0:50:150)
    ylabel(qns_roe_labels(i));
    ytickformat('%.3f')
    grid on; hold off;

    if i > 3
        xlabel('Orbit Periods')
    end
end
legend('Pre-fit','Post-fit')

%% Calculate statistics of residuals
mean_prefit = mean(prefit_resids_norm);
std_prefit = std(prefit_resids_norm);
mean_postfit = mean(postfit_resids_norm);
std_postfit = std(postfit_resids_norm);




%% Combine statistics into a table
% Rows: QNS ROEs [dsma, dlambda, dex, dey, dix, diy]
STATS_CTRL = table(ctrl_err_true_mean',ctrl_err_true_std',...
    'VariableNames',["Control Error True Mean","Control Error True Std"])
STATS_NAV = table(nav_err_mean_formal',nav_err_std_formal',nav_err_true_mean',nav_err_true_std',mean_prefit',std_prefit',mean_postfit',std_postfit',...
    'VariableNames',["Nav Error Formal Mean","Nav Error Formal Std","Nav Error True Mean", "Nav Error True Std", "Prefit Mean Residuals", "Prefit Std Residuals", "Postfit Mean Residuals", "Postfit Std Residuals"])

% Rows: RTN pos/vel [x,y,z,vx,vy,vz]
STATS_RTN = table(ctrl_err_rtn_mean',ctrl_err_rtn_std',nav_err_rtn_mean',nav_err_rtn_std',...
    'VariableNames',["Control Error True Mean","Control Error True Std","Nav Error True Mean", "Nav Error True Std"])