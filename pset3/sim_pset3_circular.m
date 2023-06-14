%% We are Close in Near-Circular Orbits
% AA 279D Problem Set 3
% Sydney Hsu and Pol Francesch
clear; clc; close all;

initial_conditions_pset3;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Parts b
% Initial HCW Orbital Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inertial relative position and velocity
rho_0 = r_init_d - r_init_c; % ECI coordinates
%norm(rho_0) / norm(r0)
rhodot_0 = v_init_d - v_init_c; % ECI coordinates / time derivative taken in ECI frame

% Relative position and velocity in the RTN frame of the chief
[rho_0_RTN,rhodot_0_RTN] = relative_motion(rho_0,rhodot_0,sma_c,ecc_c,inc_c,raan_c,qns_oe_init_c(6),mu);

state0 = [rho_0_RTN; rhodot_0_RTN];

% Check that this is valid under the HCW assumptions
% Small relative separation (rho~0.001r0)
rho_over_r0 = norm(rho_0) / norm(r_init_c);
% Equal semi-major axes
sma_diff = sma_c - oe_init_d1(1);
% Small eccentricity (e ~ 0.001)

% Orbital element difference between the deputy and the chief
oe_diff = oe_init_d1 - oe_init_c;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part c
% HCW integration constants at t0=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0 = 0;

mat1 = [sma_c, 0, 0, 0, 0, 0;
        0, sma_c, 0, 0, 0, 0;
        0, 0, sma_c, 0, 0, 0;
        0, 0, 0, sma_c*n, 0, 0;
        0, 0, 0, 0, sma_c*n, 0;
        0, 0, 0, 0, 0, sma_c*n];
mat2 = [1, sin(n*t0), cos(n*t0), 0, 0, 0;
        -3/2*n*t0, 2*cos(n*t0), -2*sin(n*t0), 1, 0, 0;
        0, 0, 0, 0, sin(n*t0), cos(n*t0);
        0, cos(n*t0), -sin(n*t0), 0, 0, 0;
        -3/2, -2*sin(n*t0), -2*cos(n*t0), 0, 0, 0;
        0, 0, 0, 0, cos(n*t0), -sin(n*t0)];
K_constants = inv(mat2)*inv(mat1)*state0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part d
% HCW solution in rectilinear coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_orbits = 50;
stepSize = T/100;
tspan = 0:stepSize:T*num_orbits;
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances
[t,state] = ode45(@(t,state) hcw(t,state,n),tspan,state0,options);

x_hcw = state(:,1); % R component
y_hcw = state(:,2); % T component
z_hcw = state(:,3); % N component
xdot_hcw = state(:,4);
ydot_hcw = state(:,5);
zdot_hcw = state(:,6);

% Position plots in the TR, NR, TN planes
figure()
% rtn_plot_pos(x_hcw, y_hcw, z_hcw);
plot_rtn_pos([x_hcw, y_hcw, z_hcw],'-',1);

% Velocity plots in the TR, NR, TN planes
figure()
% rtn_plot_vel(xdot_hcw, ydot_hcw, zdot_hcw);
plot_rtn_vel([xdot_hcw, ydot_hcw, zdot_hcw],'-',1);

% Along-track drift
figure()
orbit_periods = t/T;
plot(orbit_periods,y_hcw)
xlabel('Orbital periods')
ylabel('T (m)')




