%% EKF
% AA 279D Problem Set 8
% Sydney Hsu and Pol Francesch
clear; clc; close all;

set(0,'defaultTextInterpreter','latex');

initial_conditions_pset8;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize EKF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(1);

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
% Ground truth generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time span to validate over
orb_rev = 5;
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

% Storage variables for EKF
qns_roe_ekf  = zeros(steps+1,6);
qns_roe_ekf(1,:) = qns_roe_ekf_t;

% Store pre-fit and post-fit residuals
prefit_resids = zeros(steps+1,6);
postfit_resids = zeros(steps+1,6);
% Store the state error (covariance P)
P_diags = zeros(steps+1,6);

% STM - invariant
sma0 = oe_init_c(1); ecc0 = oe_init_c(2); inc0 = oe_init_c(3); aop0 = oe_init_c(5);
stm = stm_qns_roe_j2(stepSize,sma0,ecc0,inc0,aop0);
% stm = stm_steindorf_ms(mu, rE, J2, stepSize, sma0,ecc0,inc0,aop0);

% Looping
for j=1:steps  
    %%%% GROUND TRUTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply computed control and timestep
    tspan_local = 0:stepSize/2:stepSize;
    [~,eci_c_tspan] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_local,eci_c_t,options);
    [~,eci_d_tspan] = ode45(@(t,z) dfq_ptb(t,z,rE,mu,J2, u(j,:), eci_d_t),tspan_local,eci_d_t,options);

    eci_c_t = eci_c_tspan(end,:)'; % take the last step
    eci_d_t = eci_d_tspan(end,:)';

    [osc_sing_oe, mean_sing_oe, osc_oe,mean_oe,osc_roe,mean_roe] = rv2osc_and_mean(mu, rE, J2, 1, 0, eci_c_t, eci_d_t);

    % Save relevant variables
    eci_c(j+1,:) = eci_c_t;
    eci_d(j+1,:) = eci_d_t;
    qns_oe_c(j+1,:) = osc_oe(:,1);
    qns_oe_d(j+1,:) = osc_oe(:,2);
    qns_roe(j+1,:)  = osc_roe;

    %%%% MEASUREMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    osc_roe_meas = (sqrtm(sigma)*randn(6,1) +  (osc_roe'*sma_c))/sma_c;

    %%%% EKF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Time update
    qns_roe_minus = stm * qns_roe_ekf_t + 0; % no control
    P_minus = stm * P * stm' + Q;

    % Measurement update
    K = P_minus * H' / (H * P_minus * H' + R); % 
%     qns_roe_ekf_t = qns_roe_minus + K*(osc_roe_meas - H*qns_roe_minus);
    prefit_residual = osc_roe_meas - H*qns_roe_minus;
    qns_roe_ekf_t = qns_roe_minus + K*(prefit_residual); 
    P = (eye(6) - K*H) * P_minus * (eye(6) - K*H)' + K*R*K'; % estimate covariance

    % Save relevant variables
    qns_roe_ekf(j+1,:) = qns_roe_ekf_t;
    P_diags(j+1,:) = diag(P);
    prefit_resids(j+1,:) = prefit_residual;
    postfit_resids(j+1,:) = osc_roe_meas - H*qns_roe_ekf_t;
end

eci_c = eci_c(1:end-1,:);
eci_d = eci_d(1:end-1,:);
qns_oe_c = qns_oe_c(1:end-1,:); % osculating QNS
qns_oe_d = qns_oe_d(1:end-1,:);
qns_roe  = qns_roe(1:end-1,:);

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
plot_rtn(rho, '-', 1);

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
plot_roe_state_w_init(qns_roe_norm, roe_init1*sma_c)

% Plot EKF ROEs
figure('Name','EKF QNS ROE')
qns_roe_ekf_norm = qns_roe_ekf*sma_c;
plot_oe(orbit_periods,qns_roe_ekf_norm,qns_roe_labels);

% Plot the EKF error
figure('Name', 'ROE EKF error'); hold on;
ekf_error = (qns_roe_ekf - qns_roe) * sma_c;
qns_roe_labels = ["$a \Delta \delta a$ (m)", "$a \Delta \delta \lambda$ (m)", "$a \Delta \delta e_x$ (m)", ...
               "$a \Delta \delta e_y$ (m)", "$a \Delta \delta i_x$ (m)", "$a \Delta \delta i_y$ (m)"];
for i = 1:length(roe_init1)
    plot(orbit_periods, ekf_error(:,i))
end
legend(qns_roe_labels,'Interpreter','latex');
ylabel('EKF Error (m)');
xlabel('Orbit Periods');
grid on;
hold off;

%% Plot estimation error with standard deviations/covariance bounds, split by ROE
figure('Name', 'Estimation error with bounds'); hold on;
P_std_dev = 3*sqrt(P_diags);
plot_oe_error(orbit_periods,ekf_error,P_std_dev,qns_roe_labels)
legend('\pm \sigma','State error')

figure('Name','Steady-State')
idx = 172;
% plot_oe(orbit_periods(1:idx),P_std_dev(1:idx,:),qns_roe_labels)
plot_oe_error(orbit_periods(1:idx),ekf_error(1:idx,:),P_std_dev(1:idx,:),qns_roe_labels)
legend('\pm 3\sigma','State error')

%% Plot superimposed true and estimated states
figure('Name', 'True and Estimated States'); hold on;
plot_oe_superimposed(orbit_periods,qns_roe_norm,qns_roe_ekf_norm,qns_roe_labels);
legend('Ground truth','Estimate')

%% Plot pre-fit and post-fit residuals
figure('Name','Residuals')
plot_oe_superimposed(orbit_periods,prefit_resids.*sma_c,postfit_resids.*sma_c,qns_roe_labels)
legend('Pre-fit','Post-fit')


%% Steady-state: last orbit
mean_state_error = mean(ekf_error(end-500:end,:))
std_state_error = std(ekf_error(end-500:end,:))

