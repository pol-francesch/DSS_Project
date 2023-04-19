%% AA 279D Problem Set 2
% Sydney Hsu and Pol Francesch
clear; clc; close all;

% Constants
mu = 3.986e14;          % Earth gravitational parameter [m^3/s^2]
rE = 6378.127e3;        % Earth radius [m]
J2 = 0.00108263;        % Earth's J2 coefficient

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part a, b
% Initial Orbital Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial chief orbital parameters
sma0  = 6892.927e3;      % semi-major axis [m]
ex0   = 1e-4;            % eccentricity component in I
ey0   = 1e-4;            % eccentricity component in J
inc0  = deg2rad(97.44);  % inclination [rad]
raan0 = deg2rad(270);    % RAAN [rad]
u0    = deg2rad(0);      % argument of latitude [rad]

% Initial deputy orbital offset from chief
de = 260 / sma0;
di = 222 / sma0;
theta = deg2rad(45); % phase angle of inclination (stable at +/- 90 deg)
phase_diff = 200; % designed phase difference
dsma = 100; % displacement in semi-major axis [m]

% Calculate deputy initial orbit parameters
phi = theta - deg2rad(phase_diff); % phase angle of eccentricity
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

% Numerical integration
orb_rev = 5;
stepSize = T/100;
tspan = 0:stepSize:T*orb_rev;
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances

% Solve the EOM
[t0,z0] = ode45(@(t,z) dfq(t,z,rE,mu,0),tspan,[r0;v0],options);
[t1,z1] = ode45(@(t,z) dfq(t,z,rE,mu,0),tspan,[r1;v1],options);

% Plot to verify orbit propagation
[xE,yE,zE] = ellipsoid(0,0,0,rE,rE,rE,20);

figure(); hold on; axis equal; grid on;
set(gca,'DefaultLineLineWidth',1)
plot3(z0(:,1),z0(:,2), z0(:,3), 'm');
plot3(z1(:,1),z1(:,2), z1(:,3), 'g');
surface(xE,yE,zE, 'FaceColor','blue','FaceAlpha',.4,'EdgeColor','black','EdgeAlpha',0.5);
title('Orbit around the Earth in ECI');
xlabel('Position along I (m)');
ylabel('Position along J (m)');
zlabel('Position along K (m)')
legend('Orbit 1', 'Orbit 2','Earth');
view(3);
hold off;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part b
% Relative Numerical Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inertial relative position and velocity
rho_0 = r1 - r0; % ECI coordinates
rhodot_0 = v1 - v0; % ECI coordinates / time derivative taken in ECI frame
theta0_0 = u0; % ECI coords
[rho_0_RTN,rhodot_0_RTN] = relative_motion(rho_0,rhodot_0,sma0,ecc0,inc0,raan0,u0,mu);

% Chief pos/vel in RTN coords
[r0_0_RTN,rdot0_0_RTN] = relative_motion(r0,v0,sma0,ecc0,inc0,raan0,u0,mu);
r0_0 = norm(r0_0_RTN); % RTN coords
rdot0_0 = norm(rdot0_0_RTN); % RTN coords / RTN frame
thetadot0_0 = sqrt(mu/(sma0^3*(1-ecc0^2)^3)) * (1+ecc0*cos(u0))^2; % rotation of RTN frame in ECI (RTN coords / ECI frame)

z0 = [rho_0_RTN; rhodot_0_RTN; r0_0; theta0_0; rdot0_0; thetadot0_0]; % coordinate system: RTN

% Orbital period
T = 2*pi*sqrt(sma0^3/mu); % [sec]

% Numerical integration
orb_rev = 5;
stepSize = T/100;
tspan = 0:stepSize:T*orb_rev;
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances

% Solve the EOM
[t,z] = ode45(@(t,z) dfq_rel(t,z,mu),tspan,z0,options);

% x, y, z rel. position components in RTN coordinates
x_RTN_num = z(:,1)/sma0; xdot_RTN_num = z(:,4)/sma0;
y_RTN_num = z(:,2)/sma0; ydot_RTN_num = z(:,5)/sma0;
z_RTN_num = z(:,3)/sma0; zdot_RTN_num = z(:,6)/sma0;

orbit_periods = t / T;

rtn_plot(orbit_periods,x_RTN_num,y_RTN_num,z_RTN_num,0,1);
rtn_plot(orbit_periods,xdot_RTN_num,ydot_RTN_num,zdot_RTN_num,1,1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part c
% Absolute to Relative Analytical Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steps = length(t);
z_analytical = zeros(steps, 6);
M0 = 0;
t0 = 0;

for j=1:1:steps
    % Chief
    % Propagate mean anomaly
    M = M0 + sqrt(mu/sma0^3)*(t(j)-t0);
    M = wrapTo2Pi(M);

    % Calculate arg latitude through eccentric anomaly
    E = calc_eccentric(M, ecc0, 1e-3);
    u0 = 2 * atan(sqrt((1+ecc0)/(1-ecc0)) * tan(E/2));

    % Position & velocity in ECI
    [r0_ECI,v0_ECI] = oe2rv(mu,sma0,ex0,ey0,inc0,raan0,u0);

    % Deputy
    % Propagate mean anomaly
    M = M0 + sqrt(mu/sma1^3)*(t(j)-t0);
    M = wrapTo2Pi(M);

    % Calculate arg latitude through eccentric anomaly
    E = calc_eccentric(M, ecc1, 1e-3);
    u1 = 2 * atan(sqrt((1+ecc1)/(1-ecc1)) * tan(E/2));

    % Position & velocity in ECI
    [r1_ECI,v1_ECI] = oe2rv(mu,sma1,ex1,ey1,inc1,raan1,u1);

    % Relative
    % In ECI
    rho_ECI = r1_ECI - r0_ECI;
    rhodot_ECI = v1_ECI - v0_ECI;

    % Convert from ECI to RTN
    T_ECI2RTN = (rtn2eci_rot(inc0,raan0,u0))'; % rotation matrix from ECI to RTN
    u_dot = sqrt(mu/(sma0^3*(1-ecc0^2)^3)) * (1+ecc0*cos(u0))^2;
    w_rtn2eci = [0;0;u_dot]; % rotation of the RTN frame with respect to ECI

    rho_RTN = T_ECI2RTN*rho_ECI;
    rhodot_RTN = T_ECI2RTN*rhodot_ECI - cross(w_rtn2eci,rho_RTN);
    
    z_analytical(j,:) = [rho_RTN; rhodot_RTN]';
end

x_RTN_analytical = z_analytical(:,1)/sma0;
y_RTN_analytical = z_analytical(:,2)/sma0;
z_RTN_analytical = z_analytical(:,3)/sma0;
xdot_RTN_analytical = z_analytical(:,4)/sma0;
ydot_RTN_analytical = z_analytical(:,5)/sma0;
zdot_RTN_analytical = z_analytical(:,6)/sma0;

x_RTN = [x_RTN_num, x_RTN_analytical];
y_RTN = [y_RTN_num, y_RTN_analytical];
z_RTN = [z_RTN_num, z_RTN_analytical];
xdot_RTN = [xdot_RTN_num, xdot_RTN_analytical];
ydot_RTN = [ydot_RTN_num, ydot_RTN_analytical];
zdot_RTN = [zdot_RTN_num, zdot_RTN_analytical];

orbit_periods = t / T;

rtn_plot(orbit_periods,x_RTN,y_RTN,z_RTN,0,2);
legend('Numerical','Analytical','FontSize',12)
rtn_plot(orbit_periods, xdot_RTN, ydot_RTN, zdot_RTN, 1, 2);
legend('Numerical','Analytical','FontSize',12)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part d
% Absolute to Relative Simulation Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_error = x_RTN_num - x_RTN_analytical;
y_error = y_RTN_num - y_RTN_analytical;
z_error = z_RTN_num - z_RTN_analytical;

figure()
subplot(3,1,1);
fontSize = 11;
plot(orbit_periods,x_error);
ylabel("X-pos error (m)",'FontSize',fontSize)
title("Position error between numerical and analytical solutions",'FontSize',fontSize)
subplot(3,1,2);
plot(orbit_periods,y_error);
ylabel("Y-pos error (m)",'FontSize',fontSize)
subplot(3,1,3);
plot(orbit_periods,z_error);
ylabel("Z-pos error (m)",'FontSize',fontSize)
xlabel("Orbital periods",'FontSize',fontSize);

% Velocity errors
xdot_error = xdot_RTN_num - xdot_RTN_analytical;
ydot_error = ydot_RTN_num - ydot_RTN_analytical;
zdot_error = zdot_RTN_num - zdot_RTN_analytical;

figure()
subplot(3,1,1);
fontSize = 11;
plot(orbit_periods,xdot_error);
ylabel("X-vel error (m/s)",'FontSize',fontSize)
title("Velocity error between numerical and analytical solutions",'FontSize',fontSize)
subplot(3,1,2);
plot(orbit_periods,ydot_error);
ylabel("Y-vel error (m/s)",'FontSize',fontSize)
subplot(3,1,3);
plot(orbit_periods,zdot_error);
ylabel("Z-vel error (m/s)",'FontSize',fontSize)
xlabel("Orbital periods",'FontSize',fontSize);



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 1 - Part e / f
% Maneuver to re-establish bounded periodic motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Precalculation: delta-V required to raise semimajor axis by 100m
dsma_deltaV = -100;
v_p = sqrt(mu*(2/norm(r1)-1/sma1));
r_p = sma1*(1-ecc1);
deltaV = dsma_deltaV*mu/(2*sma1^2*v_p); % [m/s]
deltaV_vec = [0; deltaV; 0];

% Relative pos/vel in RTN coords
rho_0 = r1 - r0; % ECI coordinates
rhodot_0 = v1 - v0; % ECI coordinates / time derivative taken in ECI frame
theta0_0 = u0; % ECI coords
[rho_0_RTN,rhodot_0_RTN] = relative_motion(rho_0,rhodot_0,sma0,ecc0,inc0,raan0,u0,mu);

thetadot0_0 = sqrt(mu/(sma0^3*(1-ecc0^2)^3)) * (1+ecc0*cos(u0))^2; % RTN coords / ECI frame
w_RTNinECI = [0;0;thetadot0_0]; % angular velocity of RTN frame w.r.t. ECI expressed in the RTN frame

% Chief pos/vel in RTN coords
[r0_0_RTN,rdot0_0_RTN] = relative_motion(r0,v0,sma0,ecc0,inc0,raan0,u0,mu);
r0_0 = norm(r0_0_RTN); % RTN coords
rdot0_0 = norm(rdot0_0_RTN); % RTN coords / RTN frame

% Single orbit propagation
z0 = [rho_0_RTN; rhodot_0_RTN; r0_0; theta0_0; rdot0_0; thetadot0_0]; % coordinate system: RTN
numOrbits = 3; % number of orbits until executing maneuver
tspan_1orbit = 0:stepSize:numOrbits*T;
[t_1orbit,state_1orbit] = ode45(@(t,z) dfq_rel(t,z,mu),tspan_1orbit,z0,options);

% Apply maneuver
v1_RTN = T_ECI2RTN*(v1) - cross(w_RTNinECI,T_ECI2RTN*r1) + deltaV_vec;
v0_RTN = T_ECI2RTN*(v0) - cross(w_RTNinECI,T_ECI2RTN*r0); % Coriolis theorem
rhodot_dV_RTN = v1_RTN - v0_RTN;

rho_0_RTN = state_1orbit(end, 1:3)';

% Solve the EOM
z0_dV = [rho_0_RTN; rhodot_dV_RTN; r0_0; theta0_0; rdot0_0; thetadot0_0]; % coordinate system: RTN
[t_dV,state_dV] = ode45(@(t,z) dfq_rel(t,z,mu),tspan_1orbit,z0_dV,options);

x_1orbit = state_1orbit(1:end,1)/sma0; xdot_1orbit = state_1orbit(1:end,4)/sma0;
y_1orbit = state_1orbit(1:end,2)/sma0; ydot_1orbit = state_1orbit(1:end,5)/sma0;
z_1orbit = state_1orbit(1:end,3)/sma0; zdot_1orbit = state_1orbit(1:end,6)/sma0; 
x_dV = state_dV(:,1)/sma0;             xdot_dV = state_dV(:,4)/sma0;
y_dV = state_dV(:,2)/sma0;             ydot_dV = state_dV(:,5)/sma0;
z_dV = state_dV(:,3)/sma0;             zdot_dV = state_dV(:,6)/sma0;

x_RTN = [x_1orbit,x_dV];               xdot_RTN = [xdot_1orbit,xdot_dV];
y_RTN = [y_1orbit,y_dV];               ydot_RTN = [ydot_1orbit,ydot_dV];
z_RTN = [z_1orbit,z_dV];               zdot_RTN = [zdot_1orbit,zdot_dV];

T = 2*pi*sqrt(sma0^3/mu); % [sec]
orbit_periods_dV = [0:stepSize:numOrbits*T; ...
                    numOrbits*T:stepSize:numOrbits*2*T]'./T;

% Visualization
rtn_plot2(orbit_periods_dV, x_RTN, y_RTN, z_RTN, 0, 1);
legend('Drifting','After maneuver','Location','best','FontSize',fontSize);
rtn_plot2(orbit_periods_dV, xdot_RTN, ydot_RTN, zdot_RTN, 1, 1);
legend('Drifting','After maneuver','Location','best','FontSize',fontSize);






