%% Formation keeping simulation for mode 1 using D'Amico's thesis
% AA 279D Problem Set 5
% Sydney Hsu and Pol Francesch
clear; clc; close all;

initial_conditions_pset9;

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
e_max_man = 0.7*e_max;

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time span to validate over
orb_rev = 150;
% stepSize = T/100;
stepSize = 10.0;
% tspan = 0:60.0:3 * 86400.0;
tspan = 0:stepSize:orb_rev*T;
steps = length(tspan);

% Looping vars
roe = roe_init1;
% oe_c = oe_init_c;
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

% Impulsive control
dV = [];
t_man = [];
u_true = 0;

% Looping
for j=1:steps
    % Check if we are close to a maneuver time
    if ~isempty(t_man) && abs(t_man(1) - tspan(j)) <= stepSize
        % IP maneuvers seem to be sensitive to timing inaccuracies
        if size(dV,2) == 3
            dv_t1 = n*sma_c/4*((dsma_man - roe(1)) + norm(decc_man - roe(3:4)));
            dv_t2 = n*sma_c/4*((dsma_man - roe(1)) - norm(decc_man - roe(3:4)));
            dV(:,1) = [0; dv_t1; 0];
            dV(:,2) = [0; dv_t2; 0];
        end
        % OOP maneuvers aren't that sensitive to timing inaccuracies
        u(j,:) = dV(:,1)' / stepSize;

        t_man = t_man(2:end);
        dV    = dV(:,2:end);

        % if isempty(t_man)
        %     t_man = [0];
        %     roe(3:4)*sma_c
        % end
    end

    % Check if we are out of control window
    % Only consider this if we have no maneuvers in buffer
    % We only perform maneuvers within one/two orbit so this is good enough
    if abs(norm(roe(3:4) - roe_init1(3:4))) > e_max && isempty(t_man)
        % tspan(j)
        % Guidance - 2.4.4
        decc_man = [roe_init1(3)*cos(dphi) - roe_init1(4)*sin(dphi);
                    roe_init1(3)*sin(dphi) + roe_init1(4)*cos(dphi)];
        dinc_man = [roe_init1(5);
                    roe_init1(6) - sign(roe(5))*dinc_max];
                
        du_j2       = -12*gamma*sin(2*inc_c)*roe(5)*n*dt;
        dsma_man    = -pi/(2*n*dt - pi)*(3*e_max_man + roe(1) - 4/(3*pi)*(roe(2) - roe_init1(2) + du_j2));
        % dsma_man = 0;

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
                       u_M_oop/(2*pi)*T+T] + tspan(j) + (T - qns_oe_d(j,6) / (2*pi)*T);%mod(tspan(j), T));
        
        % Ensure OOP comes last
        while t_man_curr(3) < t_man_curr(2)
            t_man_curr(3) = t_man_curr(3) + T;
        end
        
        % Store in buffer
        dV = dV_man_curr;
        t_man = t_man_curr;
        u_M = [u_M1, u_M2, u_M_oop];
    end
    
    % Apply computed control and timestep
    tspan_local = 0:stepSize/2:stepSize;
    [~,eci_c_tspan] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan_local,eci_c_t,options);
    [~,eci_d_tspan] = ode45(@(t,z) dfq_ptb(t,z,rE,mu,J2, u(j,:), eci_d_t),tspan_local,eci_d_t,options);

    eci_c_t = eci_c_tspan(end,:)'; % take the last step
    eci_d_t = eci_d_tspan(end,:)';

    [osc_sing_oe, mean_sing_oe, osc_oe,mean_oe,osc_roe,mean_roe] = rv2osc_and_mean(mu, rE, J2, 1, 0, eci_c_t, eci_d_t);

    % Update looping vars
    % roe = mean_roe';
%     oe_c = mean_sing_oe(1,:);
    roe = (sqrtm(sigma)*randn(6,1) +  (mean_roe'*sma_c))/sma_c; % corrupt the ROEs

    % Save relevant variables
    eci_c(j+1,:) = eci_c_t;
    eci_d(j+1,:) = eci_d_t;
    qns_oe_c(j+1,:) = osc_oe(:,1);
    qns_oe_d(j+1,:) = osc_oe(:,2);
    qns_roe(j+1,:)  = osc_roe;

    u_true = u_true + (qns_oe_d(j+1,6) - qns_oe_d(j,6));
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

%% Plot initial and final states in state space
figure('Name','ROE State Space Evolution');
roe_desired = roe_init1*sma_c;
plot_roe_state_w_init_and_final(qns_roe_norm,roe_init1*sma_c,roe_desired);
%%
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
for j=1:length(roe_init1)
    plot(orbit_periods,roe_error(:,j)); hold on
end
xlabel('Orbit Periods')
ylabel('ROE error [m]')

qns_roe_labels = ["$a\delta a$", "$a\delta \lambda$", "$a\delta e_x$", ...
               "$a\delta e_y$", "$a\delta i_x$", "$a\delta i_y$"];
legend(qns_roe_labels,'Interpreter','latex');
grid on;

% Plot with control window
figure('Name', 'ROE State Space')
plot_roe_state_w_init_and_control_bounds(qns_roe_norm, roe_init1*sma_c, e_max*sma_c);
