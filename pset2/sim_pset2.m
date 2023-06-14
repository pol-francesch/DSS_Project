%% Everything is Relative
% AA 279D Problem Set 2
% Sydney Hsu and Pol Francesch
clear; clc; close all;

initial_conditions_pset2;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part b
% Relative Numerical Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inertial relative position and velocity
rho_0 = r_init_d - r_init_c; % ECI coordinates
rhodot_0 = v_init_d - v_init_c; % ECI coordinates / time derivative taken in ECI frame
theta0_0 = qns_oe_init_c(6); % ECI coords
[rho_0_RTN,rhodot_0_RTN] = relative_motion(rho_0,rhodot_0,sma_c,ecc_c,inc_c,raan_c,qns_oe_init_c(6),mu);

% Chief pos/vel in RTN coords
[r0_0_RTN,rdot0_0_RTN] = relative_motion(r_init_c,v_init_c,sma_c,ecc_c,inc_c,raan_c,oe_init_c(6),mu);
r0_0 = norm(r0_0_RTN); % RTN coords
rdot0_0 = norm(rdot0_0_RTN); % RTN coords / RTN frame
thetadot0_0 = sqrt(mu/(sma_c^3*(1-ecc_c^2)^3)) * (1+ecc_c*cos(oe_init_c(6)))^2; % rotation of RTN frame in ECI (RTN coords / ECI frame)

z0 = [rho_0_RTN; rhodot_0_RTN; r0_0; theta0_0; rdot0_0; thetadot0_0]; % coordinate system: RTN

% Numerical integration
orb_rev = 5;
stepSize = T/100;
tspan = 0:stepSize:T*orb_rev;
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances

% Solve the EOM
[t,z] = ode45(@(t,z) dfq_rel(t,z,mu),tspan,z0,options);

% Plot
orbit_periods = t / T;

figure('Name','RTN pos')
rtn_plot(orbit_periods,z(:,1),z(:,2),z(:,3),0,1);
figure('Name','RTN vel')
rtn_plot(orbit_periods,z(:,4),z(:,5),z(:,6),1,1);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part c
% Absolute to Relative Analytical Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steps = length(t);
z_analytical = zeros(steps, 6);
M0_c = oe_init_c(6);
M0_d = oe_init_d1(6);
t0 = 0;

oe_c = oe_init_c;
oe_d = oe_init_d1;

for j=1:1:steps
    % Chief
    % Propagate mean anomaly
    M = M0_c + n*(t(j)-t0);
    M = wrapTo2Pi(M);

    oe_c(6) = M;
    qns_oe_c = singular2qns(oe_c);
    u_c = qns_oe_c(6);
    uf_c = mean2true(u_c, ecc_c, 1e-8);
    f0 = mean2true(M, ecc_c, 1e-8);

    % Position & velocity in ECI
    [r0_ECI,v0_ECI] = singular_oe2rv(mu,oe_c);

    % Deputy
    % Propagate mean anomaly
    M = M0_d + n*(t(j)-t0);
    M = wrapTo2Pi(M);
    oe_d(6) = M;

    % Position & velocity in ECI
    [r1_ECI,v1_ECI] = singular_oe2rv(mu,oe_d);

    % Relative
    % In ECI
    rho_ECI = r1_ECI - r0_ECI;
    rhodot_ECI = v1_ECI - v0_ECI;

    % Convert from ECI to RTN
    T_ECI2RTN = (rtn2eci_rot(inc_c,raan_c,uf_c))'; % rotation matrix from ECI to RTN
    f_dot = sqrt(mu/(sma_c^3*(1-ecc_c^2)^3)) * (1+ecc_c*cos(uf_c))^2;
    w_rtn2eci = [0;0;f_dot]; % rotation of the RTN frame with respect to ECI

    rho_RTN = T_ECI2RTN*rho_ECI;
    rhodot_RTN = T_ECI2RTN*rhodot_ECI - cross(w_rtn2eci,rho_RTN);
    
    z_analytical(j,:) = [rho_RTN; rhodot_RTN]';
end

x_RTN = [z(:,1), z_analytical(:,1)];
y_RTN = [z(:,2), z_analytical(:,2)];
z_RTN = [z(:,3), z_analytical(:,3)];
xdot_RTN = [z(:,4), z_analytical(:,4)];
ydot_RTN = [z(:,5), z_analytical(:,5)];
zdot_RTN = [z(:,6), z_analytical(:,6)];

orbit_periods = t / T;

figure('Name','Overlay pos')
rtn_plot(orbit_periods,x_RTN,y_RTN,z_RTN,0,1);
legend('Numerical','Analytical','FontSize',12)
figure('Name','Overlay vel')
rtn_plot(orbit_periods, xdot_RTN, ydot_RTN, zdot_RTN, 1, 1);
legend('Numerical','Analytical','FontSize',12)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part d
% Absolute to Relative Simulation Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_error = z(:,1) - z_analytical(:,1);
y_error = z(:,2) - z_analytical(:,2);
z_error = z(:,3) - z_analytical(:,3);

figure('Name','RTN position error')
subplot(3,1,1);
% fontSize = 20;
plot(orbit_periods,x_error); xticks(1:1:orb_rev);
ylabel("R (m)");
% title("Position error between numerical and analytical solutions",'FontSize',fontSize)
subplot(3,1,2);
plot(orbit_periods,y_error); xticks(1:1:orb_rev);
ylabel("T (m)");
subplot(3,1,3);
plot(orbit_periods,z_error); xticks(1:1:orb_rev);
ylabel("N (m)");
xlabel("Orbital periods");

% Velocity errors
xdot_error = z(:,4) - z_analytical(:,4);
ydot_error = z(:,5) - z_analytical(:,5);
zdot_error = z(:,6) - z_analytical(:,6);

figure('Name','RTN velocity error')
subplot(3,1,1);
% fontSize = 11;
plot(orbit_periods,xdot_error); xticks(1:1:orb_rev);
ylabel("R (m/s)");
% title("Velocity error between numerical and analytical solutions",'FontSize',fontSize)
subplot(3,1,2);
plot(orbit_periods,ydot_error); xticks(1:1:orb_rev);
ylabel("T (m/s)");
subplot(3,1,3);
plot(orbit_periods,zdot_error); xticks(1:1:orb_rev);
ylabel("N (m/s)");
xlabel("Orbital periods");


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part e / f
% Maneuver to re-establish bounded periodic motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Precalculation: delta-V required to raise semimajor axis by 100m
dsma_deltaV = -100; % make sure this matches the initial conditions
v_p = sqrt(mu*(2/norm(r_init_d)-1/oe_init_d1(1)));
r_p = oe_init_d1(1)*(1-oe_init_d1(2));
deltaV = dsma_deltaV*mu/(2*oe_init_d1(1)^2*v_p); % [m/s]
deltaV_vec = [0; deltaV; 0];

% Inertial relative position and velocity
rho_0 = r_init_d - r_init_c; % ECI coordinates
rhodot_0 = v_init_d - v_init_c; % ECI coordinates / time derivative taken in ECI frame
theta0_0 = qns_oe_init_c(6); % ECI coords
[rho_0_RTN,rhodot_0_RTN] = relative_motion(rho_0,rhodot_0,sma_c,ecc_c,inc_c,raan_c,qns_oe_init_c(6),mu);

thetadot0_0 = sqrt(mu/(sma_c^3*(1-ecc_c^2)^3)) * (1+ecc_c*cos(oe_init_c(6)))^2;
w_RTNinECI = [0;0;thetadot0_0]; % angular velocity of RTN frame w.r.t. ECI expressed in the RTN frame

% Chief pos/vel in RTN coords
[r0_0_RTN,rdot0_0_RTN] = relative_motion(r_init_c,v_init_c,sma_c,ecc_c,inc_c,raan_c,oe_init_c(6),mu);
r0_0 = norm(r0_0_RTN); % RTN coords
rdot0_0 = norm(rdot0_0_RTN); % RTN coords / RTN frame

% Three orbit propagation
z0 = [rho_0_RTN; rhodot_0_RTN; r0_0; theta0_0; rdot0_0; thetadot0_0]; % coordinate system: RTN
numOrbits = 3; % number of orbits until executing maneuver
tspan_1orbit = 0:stepSize:numOrbits*T;
[t_1orbit,state_1orbit] = ode45(@(t,z) dfq_rel(t,z,mu),tspan_1orbit,z0,options);

% Apply maneuver
v1_RTN = T_ECI2RTN*(v_init_d) - cross(w_RTNinECI,T_ECI2RTN*r_init_d) + deltaV_vec;
v0_RTN = T_ECI2RTN*(v_init_c) - cross(w_RTNinECI,T_ECI2RTN*r_init_c); % Coriolis theorem
rhodot_dV_RTN = v1_RTN - v0_RTN;

rho_0_RTN = state_1orbit(end, 1:3)';

% Solve the EOM
z0_dV = [rho_0_RTN; rhodot_dV_RTN; r0_0; theta0_0; rdot0_0; thetadot0_0]; % coordinate system: RTN
[t_dV,state_dV] = ode45(@(t,z) dfq_rel(t,z,mu),tspan_1orbit,z0_dV,options);

% Plot
x_RTN = [state_1orbit(1:end,1),state_dV(:,1)];               
y_RTN = [state_1orbit(1:end,2),state_dV(:,2)];               
z_RTN = [state_1orbit(1:end,3),state_dV(:,3)];             

xdot_RTN = [state_1orbit(1:end,4),state_dV(:,4)];
ydot_RTN = [state_1orbit(1:end,5),state_dV(:,5)];
zdot_RTN = [state_1orbit(1:end,6),state_dV(:,6)];

orbit_periods_dV = [0:stepSize:numOrbits*T; ...
                    numOrbits*T:stepSize:numOrbits*2*T]'./T;

figure('Name','Maneuver RTN pos')
rtn_plot2(orbit_periods_dV, x_RTN, y_RTN, z_RTN, 0, 1);
legend('Drifting','After maneuver','Location','best');

figure('Name','Maneuver RTN vel')
rtn_plot2(orbit_periods_dV, xdot_RTN, ydot_RTN, zdot_RTN, 1, 1);
legend('Drifting','After maneuver','Location','best');

