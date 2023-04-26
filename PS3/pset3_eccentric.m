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
% Problem 2 - Part a
% Initial TH Orbital Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rp = 6892.927e3;      % periapsis [m]
ecc0 = 0.1;         % eccentricity

% Initial chief orbital parameters
sma0  = rp / (1-ecc0);   % semi-major axis [m]
inc0  = deg2rad(97.44);  % inclination [rad]
raan0 = deg2rad(270);    % RAAN [rad]
aop0  = deg2rad(0);      % argument of periapsis [rad]
ta0    = deg2rad(0);     % true anomaly [rad]

% Initial deputy orbital offset from chief
de = 260 / sma0;
di = 222 / sma0;
theta = deg2rad(45); % phase angle of inclination (stable at +/- 90 deg)
phase_diff = 200; % designed phase difference
dsma = 0; %-10000; % semi-major axis difference [m]

% Calculate deputy initial orbit parameters
phi = theta - deg2rad(phase_diff); % phase angle of eccentricity
dinc = di*cos(theta);
draan = di*sin(theta) / sin(inc0);

sma1 = sma0 + dsma; ecc1 = ecc0 + de; inc1 = inc0 + dinc; 
raan1 = raan0 + draan; aop1 = aop0; ta1 = ta0;

% Initial position and velocity in inertial frame
[r0,v0] = classic_oe2rv(mu,sma0,ecc0,inc0,raan0,aop0,ta0);
[r1,v1] = classic_oe2rv(mu,sma1,ecc1,inc1,raan1,aop1,ta1);

% Orbital period
T = 2*pi*sqrt(sma0^3/mu); % [sec]

% Inertial relative position and velocity
rho_0 = r1 - r0; % ECI coordinates
% norm(rho_0) / rp
rhodot_0 = v1 - v0; % ECI coordinates / time derivative taken in ECI frame

% Relative position and velocity in the RTN frame of the chief
% Coriolis theorem
thetadot0_0 = sqrt(mu/(sma0^3*(1-ecc0^2)^3)) * (1+ecc0*cos(ta0))^2; % RTN coords / ECI frame
w_RTNinECI = [0;0;thetadot0_0]; % angular velocity of RTN frame w.r.t. ECI expressed in the RTN frame
T_ECI2RTN = (rtn2eci_rot(inc0,raan0,aop0+ta0))'; % rotation matrix from ECI to RTN
rho_0_RTN = T_ECI2RTN*rho_0;
rhodot_0_RTN = T_ECI2RTN*rhodot_0 - cross(w_RTNinECI,rho_0_RTN); % RTN coords / RTN frame

% State vector at t=0 [x; y; z; xdot; ydot; zdot] in RTN coords / RTN frame
state0 = [rho_0_RTN; rhodot_0_RTN];

% Singular ROE's
roe = [(sma1-sma0)/sma0; 0; ecc1 - ecc0; aop1 - aop0; inc1 - inc0; raan1 - raan0];


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part b
% YH integration constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stm = dimensional_ya_stm(0, mu, sma0, ecc0, ta0);
c = inv(stm)*state0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Parts c, d
% Yamanaka-Ankersen analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_orbits = 15;
stepSize = 2*pi/1000;
tspan = 0:T/1000:num_orbits*T;
anglespan = 0:stepSize:2*pi*num_orbits;
t0 = 0;
n = sqrt(mu/sma0^3);

% True anomaly as independent variable
states_c = zeros(length(tspan),6);
h = norm(cross(r0,v0));

for j=1:length(tspan)
% %     f = anglespan(j); % true anomaly
% %     u = f+aop0; % argument of latitude
% %     t = u/n; % time
    % True anomaly
    t = tspan(j);
    M = n*(t - t0);
% %     E = calc_eccentric(M, ecc0, 1e-5);
    if abs(abs(calc_eccentric(M, ecc0, 1e-5)) - abs(M)) < 1
        E = calc_eccentric(M, ecc0, 1e-5);
    end
    f = 2 * atan(sqrt((1+ecc0)/(1-ecc0)) * tan(E/2));

    % Relative position and velocity
    states_c(j,:) = dimensional_ya_stm(t, mu, sma0, ecc0, f) * c / sma0;
end

% Unbounded motion
orbit_periods = tspan / T;
figure();
subplot(2,1,1);
plot(orbit_periods, states_c(:,1));
xlabel("Orbital Periods"); ylabel("$x / a_0$"); grid on;

subplot(2,1,2);
plot(orbit_periods, states_c(:,2));
xlabel("Orbital Periods"); ylabel("$y / a_0$"); grid on;

% Position
figure();
subplot(2,2,1)
plot(states_c(:,2),states_c(:,1))
title('TR'); xlabel('$y / a_0$'); ylabel('$x / a_0$'); axis equal; grid on
subplot(2,2,2)
plot(states_c(:,3),states_c(:,1))
title('NR'); xlabel('$z / a_0$'); ylabel('$x / a_0$'); axis equal; grid on
subplot(2,2,3)
plot(states_c(:,2),states_c(:,3))
title('TN'); xlabel('$y /a_0$'); ylabel('$z /a_0$'); axis equal; grid on
subplot(2,2,4)
plot3(states_c(:,1),states_c(:,2),states_c(:,3)); 
xlabel('$x /a_0$'); ylabel('$y /a_0$'); zlabel('$z /a_0$'); grid on

% Velocity
figure();
subplot(2,2,1)
plot(states_c(:,5),states_c(:,4))
title('TR'); xlabel('$\dot{y} / a_0$'); ylabel('$\dot{x} / a_0$'); axis equal; grid on
subplot(2,2,2)
plot(states_c(:,6),states_c(:,4))
title('NR'); xlabel('$\dot{z} / a_0$'); ylabel('$\dot{x} / a_0$'); axis equal; grid on
subplot(2,2,3)
plot(states_c(:,5),states_c(:,6))
title('TN'); xlabel('$\dot{y} /a_0$'); ylabel('$\dot{z} /a_0$'); axis equal; grid on
subplot(2,2,4)
plot3(states_c(:,4),states_c(:,5),states_c(:,6)); 
xlabel('$\dot{x} /a_0$'); ylabel('$\dot{y} /a_0$'); zlabel('$\dot{z} /a_0$'); grid on

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part e
% Quasi-nonsingular relative orbital elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eccentricity vector components
ex0 = ecc0*cos(aop0); ey0 = ecc0*sin(aop0);
ex1 = ecc1*cos(aop1); ey1 = ecc1*sin(aop1);
M0 = ta2ma(aop0,ecc0);
M1 = ta2ma(aop1,ecc1);

% Quasi-nonsingular relative orbital elements between chief and deputy
da = (sma1-sma0)/sma0;
dex = ex1-ex0; 
dey = ey1-ey0;
dlambda = (M1-M0)+(aop1-aop0)+(raan1-raan0)*cos(inc0);
dix = inc1-inc0; 
diy = (raan1-raan0)*sin(inc0);

qns_roe = [da; dlambda; dex; dey; dix; diy];


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part f
% ROE geometric linear mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% True anomaly as independent variable
eta = sqrt(1-ecc0^2);
state_map = zeros(6, length(tspan));

for j=1:length(tspan)
    % True anomaly
    t = tspan(j);
    M = n*(t - t0);
% %     E = calc_eccentric(M, ecc0, 1e-5);
    if abs(abs(calc_eccentric(M, ecc0, 1e-5)) - abs(M)) < 1
        E = calc_eccentric(M, ecc0, 1e-5);
    end
    f = 2 * atan(sqrt((1+ecc0)/(1-ecc0)) * tan(E/2));
    u = f + aop0;
% %     f = anglespan(j); % true anomaly
% %     u = f+aop0; % argument of latitude
% %     t = u/n; % time
    
    % Auxiliary
    k = 1+ex0*cos(u) + ey0*sin(u);
    k_prime = -ex0*sin(u) + ey0*cos(u);
    
    % Coefficients
    bx1 = 1/k + 3/2*k_prime*n/eta^3*t;
    bx2 = -k_prime/eta^3;
    bx3 = 1/eta^3*(ex0*(k-1)/(1+eta)-cos(u));
    bx4 = 1/eta^3*(ey0*(k-1)/(1+eta)-sin(u));
    bx6 = k_prime/eta^3*cot(inc0);
    by1 = -3/2*k*n/eta^3*t;
    by2 = k/eta^3;
    by3 = 1/eta^2*((1+1/k)*sin(u)+ey0/k+k/eta*(ey0/(1+eta)));
    by4 = -1/eta^2*((1+1/k)*cos(u)+ex0/k+k/eta*(ex0/(1+eta)));
    by6 = (1/k-k/eta^3)*cot(inc0);
    bz5 = 1/k*sin(u);
    bz6 = -1/k*cos(u);
    bxd1 = k_prime/2+3/2*k^2*(1-k)*n/eta^3*t;
    bxd2 = k^2/eta^3*(k-1);
    bxd3 = k^2/eta^3*(eta*sin(u)+ey0*(k-1)/(1+eta));
    bxd4 = -k^2/eta^3*(eta*cos(u)+ex0*(k-1)/(1+eta));
    bxd6 = -k^2/eta^3*(k-1)*cot(inc0);
    byd1 = -3/2*k*(1+k*k_prime*n/eta^3*t);
    byd2 = k^2/eta^3*k_prime;
    byd3 = (1+k^2/eta^3)*cos(u)+ex0*k/eta^2*(1+k/eta*(1-k)/(1+eta));
    byd4 = (1+k^2/eta^3)*sin(u)+ey0*k/eta^2*(1+k/eta*(1-k)/(1+eta));
    byd6 = -(1+k^2/eta^3)*k_prime*cot(inc0);
    bzd5 = cos(u)+ex0;
    bzd6 = sin(u)+ey0;

    % STM
    b_matrix = [bx1, bx2, bx3, bx4, 0, bx6;
                by1, by2, by3, by4, 0, by6;
                0, 0, 0, 0, bz5, bz6;
                bxd1, bxd2, bxd3, bxd4, 0, bxd6;
                byd1, byd2, byd3, byd4, 0, byd6;
                0, 0, 0, 0, bzd5, bzd6];
    state_map(:,j) = [sma0*eta^2*eye(3), zeros(3,3); zeros(3,3), sma0*n/eta*eye(3)] * b_matrix * qns_roe;
end

x_bar_map = state_map(1,:) ./ sma0; xdot_bar_map = state_map(4,:) ./ sma0;
y_bar_map = state_map(2,:) ./ sma0; ydot_bar_map = state_map(5,:) ./ sma0;
z_bar_map = state_map(3,:) ./ sma0; zdot_bar_map = state_map(6,:) ./ sma0;

% Position plots overlaid with YA solution
figure()
subplot(2,2,1); hold on
plot(states_c(:,2),states_c(:,1),'linewidth',2); hold on; plot(y_bar_map,x_bar_map,'--');
title('TR'); xlabel('$y/a_0$'); ylabel('$x/a_0$'); axis equal; grid on
subplot(2,2,2)
plot(states_c(:,3),states_c(:,1),'linewidth',2); hold on; plot(z_bar_map,x_bar_map,'--');
title('NR'); xlabel('$z/a_0$'); ylabel('$x/a_0$'); axis equal; grid on
subplot(2,2,3)
plot(states_c(:,2),states_c(:,3),'linewidth',2); hold on; plot(y_bar_map,z_bar_map,'--');
title('TN'); xlabel('$y/a_0$'); ylabel('$z/a_0$'); axis equal; grid on
subplot(2,2,4); hold on; view(3);
plot3(states_c(:,1),states_c(:,2),states_c(:,3),'--'); xlabel('$x/a_0$'); ylabel('$y/a_0$'); zlabel('$z/a_0$'); hold on;
plot3(x_bar_map,y_bar_map,z_bar_map, '--','linewidth',2);
grid on

% Velocity plots overlaid with YA solution
figure()
subplot(2,2,1); hold on
plot(states_c(:,5),states_c(:,4),'linewidth',2); hold on; plot(ydot_bar_map,xdot_bar_map,'--');
title('TR'); xlabel('$\dot{y}/a_0$'); ylabel('$\dot{x}/a_0$'); axis equal; grid on
subplot(2,2,2)
plot(states_c(:,6),states_c(:,4),'linewidth',2); hold on; plot(zdot_bar_map,xdot_bar_map,'--');
title('NR'); xlabel('$\dot{z}/a_0$'); ylabel('$\dot{x}/a_0$'); axis equal; grid on
subplot(2,2,3)
plot(states_c(:,5),states_c(:,6),'linewidth',2); hold on; plot(ydot_bar_map,zdot_bar_map,'--');
title('TN'); xlabel('$\dot{y}/a_0$'); ylabel('$\dot{z}/a_0$'); axis equal; grid on
subplot(2,2,4); hold on; view(3);
plot3(states_c(:,4),states_c(:,5),states_c(:,6),'linewidth',2); xlabel('$\dot{x}/a_0$'); ylabel('$\dot{y}/a_0$'); zlabel('$\dot{z}/a_0$'); hold on;
plot3(xdot_bar_map,ydot_bar_map,zdot_bar_map, '--','linewidth',2);
grid on

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part g
%  Compare integration constants to ROE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
icu_2_roe = [1, 0, 0, 0, 0, 0;
             0, -ex0*(eta + 1/(1+eta)), ey0*(eta + 1/(1+eta)), 1, 0, 0;
             0, ex0*ey0, ex0^2-1, -ey0, 0, -ey0*cot(inc0);
             0, ey0^2-1, ex0*ey0, ex0, 0, ex0*cot(inc0);
             0, 0, 0, 0, 1, 0;
             0, 0, 0, 0, 0, -1];
roe_from_ic = icu_2_roe * c;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part h
%  True relative motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inertial relative position and velocity
rho_0 = r1 - r0; % ECI coordinates
rhodot_0 = v1 - v0; % ECI coordinates / time derivative taken in ECI frame
u0 = 0;
theta0_0 = u0; % ECI coords
[rho_0_RTN,rhodot_0_RTN] = relative_motion(rho_0,rhodot_0,sma0,ecc0,inc0,raan0,u0,mu);

% Chief pos/vel in RTN coords
[r0_0_RTN,rdot0_0_RTN] = relative_motion(r0,v0,sma0,ecc0,inc0,raan0,u0,mu);
r0_0 = norm(r0_0_RTN); % RTN coords
rdot0_0 = norm(rdot0_0_RTN); % RTN coords / RTN frame
thetadot0_0 = sqrt(mu/(sma0^3*(1-ecc0^2)^3)) * (1+ecc0*cos(u0))^2; % rotation of RTN frame in ECI (RTN coords / ECI frame)

z0 = [rho_0_RTN; rhodot_0_RTN; r0_0; theta0_0; rdot0_0; thetadot0_0]; % coordinate system: RTN

% Numerical integration
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances

% Solve the EOM
%%tspan = linspace(0,num_orbits*T,length(tspan));
[t,state_num] = ode45(@(t,z) dfq_rel(t,z,mu),tspan,z0,options);

x_RTN_num = state_num(:,1)/sma0; xdot_RTN_num = state_num(:,4)/sma0;
y_RTN_num = state_num(:,2)/sma0; ydot_RTN_num = state_num(:,5)/sma0;
z_RTN_num = state_num(:,3)/sma0; zdot_RTN_num = state_num(:,6)/sma0;

% Position plots overlaid with YA solution & geometric mapping
figure()
subplot(2,2,1); hold on
plot(states_c(:,2),states_c(:,1),'linewidth',3); hold on; 
plot(y_bar_map,x_bar_map,'--','linewidth',2); hold on; plot(y_RTN_num, x_RTN_num,':');
title('TR'); xlabel('$y/a_0$'); ylabel('$x/a_0$'); axis equal; grid on
subplot(2,2,2)
plot(states_c(:,3),states_c(:,1),'linewidth',3); hold on; 
plot(z_bar_map,x_bar_map,'--','linewidth',2); plot(z_RTN_num, x_RTN_num,':');
title('NR'); xlabel('$z/a_0$'); ylabel('$x/a_0$'); axis equal; grid on
subplot(2,2,3)
plot(states_c(:,2),states_c(:,3),'linewidth',3); hold on; 
plot(y_bar_map,z_bar_map,'--','linewidth',2); plot(y_RTN_num, z_RTN_num,':');
title('TN'); xlabel('$y/a_0$'); ylabel('$z/a_0$'); axis equal; grid on
subplot(2,2,4); hold on; view(3);
plot3(states_c(:,1),states_c(:,2),states_c(:,3),'linewidth',1); xlabel('$x/a_0$'); ylabel('$y/a_0$'); zlabel('$z/a_0$'); hold on;
plot3(x_bar_map,y_bar_map,z_bar_map, '--','linewidth',3); plot3(x_RTN_num, y_RTN_num, z_RTN_num,':','linewidth',2);
grid on

% Velocity plots overlaid with YA solution & geometric mapping
figure()
subplot(2,2,1); hold on
plot(states_c(:,5),states_c(:,4),'linewidth',3); hold on; 
plot(ydot_bar_map,xdot_bar_map,'--','linewidth',2); plot(ydot_RTN_num, xdot_RTN_num,':');
title('TR'); xlabel('$\dot{y}/a_0$'); ylabel('$\dot{x}/a_0$'); axis equal; grid on
subplot(2,2,2)
plot(states_c(:,6),states_c(:,4),'linewidth',3); hold on; 
plot(zdot_bar_map,xdot_bar_map,'--','linewidth',2); plot(zdot_RTN_num, xdot_RTN_num,':');
title('NR'); xlabel('$\dot{z}/a_0$'); ylabel('$\dot{x}/a_0$'); axis equal; grid on
subplot(2,2,3)
plot(states_c(:,5),states_c(:,6),'linewidth',3); hold on; 
plot(ydot_bar_map,zdot_bar_map,'--','linewidth',2); plot(ydot_RTN_num, zdot_RTN_num,':');
title('TN'); xlabel('$\dot{y}/a_0$'); ylabel('$\dot{z}/a_0$'); axis equal; grid on
subplot(2,2,4); hold on; view(3);
plot3(states_c(:,4),states_c(:,5),states_c(:,6)); xlabel('$\dot{x}/a_0$'); ylabel('$\dot{y}/a_0$'); zlabel('$\dot{z}/a_0$'); hold on;
plot3(xdot_bar_map,ydot_bar_map,zdot_bar_map, '--','linewidth',2); plot3(xdot_RTN_num, ydot_RTN_num, zdot_RTN_num,':');
grid on

% Error from true numerical propagation
ya_error = states_c - state_num(:,1:6) ./ sma0;         % YA solution error
gm_error = (state_map' - state_num(:,1:6)) ./sma0;      % Geometric mapping error

figure();
subplot(3,2,1); hold on; grid on;
plot(orbit_periods, ya_error(:,1)); 
plot(orbit_periods, gm_error(:,1));
ylabel('$x-error / a_0$'); hold off;

subplot(3,2,3); hold on; grid on;
plot(orbit_periods, ya_error(:,2)); 
plot(orbit_periods, gm_error(:,2));
ylabel('$y-error / a_0$'); hold off;

subplot(3,2,5); hold on; grid on;
plot(orbit_periods, ya_error(:,3)); 
plot(orbit_periods, gm_error(:,3));
ylabel('$z-error / a_0$'); xlabel("Orbital Periods"); hold off;

subplot(3,2,2); hold on; grid on;
plot(orbit_periods, ya_error(:,4)); 
plot(orbit_periods, gm_error(:,4));
ylabel('$\dot{x}-error / a_0$'); hold off;

subplot(3,2,4); hold on; grid on;
plot(orbit_periods, ya_error(:,5)); 
plot(orbit_periods, gm_error(:,5));
ylabel('$\dot{y}-error / a_0$'); hold off;

subplot(3,2,6); hold on; grid on;
plot(orbit_periods, ya_error(:,6)); 
plot(orbit_periods, gm_error(:,6));
ylabel('$\dot{z}-error / a_0$'); xlabel("Orbital Periods"); hold off;


