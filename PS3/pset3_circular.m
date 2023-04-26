%% AA 279D Problem Set 3
% Sydney Hsu and Pol Francesch
clear; clc; close all;

% Constants
mu = 3.986e14;          % Earth gravitational parameter [m^3/s^2]
rE = 6378.127e3;        % Earth radius [m]
J2 = 0.00108263;        % Earth's J2 coefficient
set(0,'defaultTextInterpreter','latex');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Parts a, b
% Initial HCW Orbital Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial chief orbital parameters
sma0  = 6892.927e3;      % semi-major axis [m]
ex0   = 1e-4;            % eccentricity component in X
ey0   = 1e-4;            % eccentricity component in Y
inc0  = deg2rad(97.44);  % inclination [rad]
raan0 = deg2rad(270);    % RAAN [rad]
u0    = deg2rad(0);      % argument of latitude [rad]

% Initial deputy orbital offset from chief
de = 260 / sma0;
di = 222 / sma0;
theta = deg2rad(45); % phase angle of inclination (stable at +/- 90 deg)
%phi = deg2rad(0); % Design this parameter instead?
phase_diff = 200; % designed phase difference
dsma = 0; % semi-major axis difference [m]

% Calculate deputy initial orbit parameters
% TODO: Check whether we want to design phi or theta
phi = theta - deg2rad(phase_diff); % phase angle of eccentricity
%theta = phi + deg2rad(phase_diff);
dex = de*cos(phi);  dey = de*sin(phi);
dinc = di*cos(theta);
draan = di*sin(theta) / sin(inc0);

sma1 = sma0 + dsma; ex1 = ex0 + dex; ey1 = ey0 + dey;
inc1 = inc0 + dinc; raan1 = raan0 + draan; u1 = u0;

ecc0 = sqrt(ex0^2 + ey0^2); % eccentricity magnitude
ecc1 = sqrt(ex1^2 + ey1^2);

% Initial position and velocity in inertial frame
[r0,v0] = oe2rv(mu,sma0,ex0,ey0,inc0,raan0,u0);
[r1,v1] = oe2rv(mu,sma1,ex1,ey1,inc1,raan1,u1);

% Orbital period
T = 2*pi*sqrt(sma0^3/mu); % [sec]

% Inertial relative position and velocity
rho_0 = r1 - r0; % ECI coordinates
%norm(rho_0) / norm(r0)
rhodot_0 = v1 - v0; % ECI coordinates / time derivative taken in ECI frame

% Relative position and velocity in the RTN frame of the chief
[rho_0_RTN,rhodot_0_RTN] = relative_motion(rho_0,rhodot_0,sma0,ecc0,inc0,raan0,u0,mu);
% State vector at t=0 [x; y; z; xdot; ydot; zdot] in RTN coords / RTN frame
state0 = [rho_0_RTN(1); rho_0_RTN(2); rho_0_RTN(3); rhodot_0_RTN(1); rhodot_0_RTN(2); rhodot_0_RTN(3)];

% Quasi-nonsingular relative orbital elements between chief and deputy
%roe = [(sma1-sma0)/sma0; (u1-u0)+(raan1-raan0)*cos(inc0); ex1-ex0; ey1-ey0; inc1-inc0; (raan1-raan0)*sin(inc0)];

% Orbit element differences
delta_sma = sma1-sma0;
delta_ex = ex1-ex0;
delta_ey = ey1-ey0;
delta_inc = inc1-inc0;
delta_raan = raan1-raan0;
delta_u = u1-u0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part c
% HCW integration constants at t0=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = sqrt(mu/sma0^3);
t0 = 0;

mat1 = [sma0, 0, 0, 0, 0, 0;
        0, sma0, 0, 0, 0, 0;
        0, 0, sma0, 0, 0, 0;
        0, 0, 0, sma0*n, 0, 0;
        0, 0, 0, 0, sma0*n, 0;
        0, 0, 0, 0, 0, sma0*n];
mat2 = [1, sin(n*t0), cos(n*t0), 0, 0, 0;
        -3/2*n*t0, 2*cos(n*t0), -2*sin(n*t0), 1, 0, 0;
        0, 0, 0, 0, sin(n*t0), cos(n*t0);
        0, cos(n*t0), -sin(n*t0), 0, 0, 0;
        -3/2, -2*sin(n*t0), -2*cos(n*t0), 0, 0, 0;
        0, 0, 0, 0, cos(n*t0), -sin(n*t0)];
K_constants = inv(mat2)*inv(mat1)*state0;

% Alternative formulations
% K_constants = sma0*[dsma;
%                    (u1-u0) + (raan1-raan0)*cos(inc0);
%                    dex;
%                    dey;
%                    inc1-inc0;
%                    (raan1-raan0)*sin(inc0)];

% Initial pos/vel at t0 = 0
% x0 = rho_0_RTN(1);
% y0 = rho_0_RTN(2);
% z0 = rho_0_RTN(3);
% xdot0 = rhodot_0_RTN(1);
% ydot0 = rhodot_0_RTN(2);
% zdot0 = rhodot_0_RTN(3);
% K_constants2 = [4*x0 + 2*ydot0/n;
%                 y0-2*xdot0/n;
%                 3*x0+2*ydot0/n;
%                 -xdot0/n;
%                 zdot0/n;
%                 z0];

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
subplot(2,2,1)
plot(y_hcw,x_hcw)
title('TR'); xlabel('y (m)'); ylabel('x (m)'); axis equal; grid on
subplot(2,2,2)
plot(z_hcw,x_hcw)
title('NR'); xlabel('z (m)'); ylabel('x (m)'); axis equal; grid on
subplot(2,2,3)
plot(y_hcw,z_hcw)
title('TN'); xlabel('y (m)'); ylabel('z (m)'); axis equal; grid on
subplot(2,2,4)
plot3(x_hcw,y_hcw,z_hcw); xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)'); grid on

% Velocity plots in the TR, NR, TN planes
figure()
subplot(2,2,1)
plot(ydot_hcw,xdot_hcw)
title('TR'); xlabel('$\dot{y}$ (m/s)'); ylabel('$\dot{x}$ (m/s)'); axis equal; grid on
subplot(2,2,2)
plot(zdot_hcw,xdot_hcw)
title('NR'); xlabel('$\dot{z}$ (m/s)'); ylabel('$\dot{x}$ (m/s)'); axis equal; grid on
subplot(2,2,3)
plot(ydot_hcw,zdot_hcw)
title('TN'); xlabel('$\dot{y}$ (m/s)'); ylabel('$\dot{z}$ (m/s)'); axis equal; grid on
subplot(2,2,4)
plot3(xdot_hcw,ydot_hcw,zdot_hcw); xlabel('$\dot{x}$ (m/s)'); ylabel('$\dot{y}$ (m/s)'); zlabel('$\dot{z}$ (m/s)'); grid on

% Along-track drift
figure()
orbit_periods = t/T;
plot(orbit_periods,y_hcw)
xlabel('Orbit periods')
ylabel('y (m)')



