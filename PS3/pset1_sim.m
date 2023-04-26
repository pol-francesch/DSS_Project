%% AA 279D Problem Set 1
% Sydney Hsu and Pol Francesch
clear; clc; close all;

% Constants
mu = 3.986e14;          % Earth gravitational parameter [m^3/s^2]
rE = 6378.127e3;        % Earth radius [m]
J2 = 0.00108263;        % Earth's J2 coefficient

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part a, b, c
% Numerical Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% New initial orbital parameters
sma0  = 6892.927e3;      % semi-major axis [m]
ex0   = 1e-4;            % eccentricity component in I
ey0   = 1e-4;            % eccentricity component in J
inc0  = deg2rad(97.44);  % inclination [rad]
raan0 = deg2rad(270);    % RAAN [rad]
u0    = deg2rad(0);      % argument of latitude [rad]

sma = sma0; ex = ex0; ey = ey0;
inc = inc0; raan = raan0; u = u0;

ecc0 = sqrt(ex^2 + ey^2);
ecc = ecc0;

% Initial position and velocity in inertial frame
[r0,v0] = oe2rv(mu,sma,ex,ey,inc,raan,u);

% Orbital period
T = 2*pi*sqrt(sma^3/mu); % [sec]

% Numerical integration
orb_rev = 5;
stepSize = T/100;
tspan = 0:stepSize:T*orb_rev;
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances
z0 = [r0;v0];

% Solve the EOM
[t,z]       = ode45(@(t,z) dfq(t,z,rE,mu,0),tspan,z0,options);
[t_j2,z_j2] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan,z0,options);

% Plot
[xE,yE,zE] = ellipsoid(0,0,0,rE,rE,rE,20);

figure(); hold on; axis equal; grid on;
set(gca,'DefaultLineLineWidth',1)
plot3(z(:,1),z(:,2), z(:,3), 'm');
surface(xE,yE,zE, 'FaceColor','blue','FaceAlpha',.4,'EdgeColor','black','EdgeAlpha',0.5);
title('Orbit around the Earth in ECI');
xlabel('Position along I (m)');
ylabel('Position along J (m)');
zlabel('Position along K (m)')
legend('Orbit','Earth');
view(3);
hold off;

figure(); hold on; axis equal; grid on;
set(gca,'DefaultLineLineWidth',1)
plot3(z_j2(:,1),z_j2(:,2), z_j2(:,3), 'm');
surface(xE,yE,zE, 'FaceColor','blue','FaceAlpha',.4,'EdgeColor','black','EdgeAlpha',0.5);
title('Perturbed Orbit around the Earth in ECI');
xlabel('Position along I (m)');
ylabel('Position along J (m)');
zlabel('Position along K (m)')
legend('Orbit','Earth');
view(3);
hold off;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part d
% Comparing analytical propagation to numerical
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analytical propagation in ECI
steps = length(t);

z_analytic = zeros(steps,6); % position in the ECI frame
M0 = 0;
t0 = 0;

pos_err = zeros(steps,3);
vel_err = zeros(steps,3);

for j = 1:steps
    % Propagate mean anomaly
    M = M0 + sqrt(mu/sma^3)*(t(j)-t0);
    M = wrapTo2Pi(M);

    % Calculate true anomaly through eccentric anomaly
    E = calc_eccentric(M, ecc, 1e-3);
    u = 2 * atan(sqrt((1+ecc)/(1-ecc)) * tan(E/2));
    
    % Position & velocity in ECI
    [r_analytic_ECI,v_analytic_ECI] = oe2rv(mu,sma,ex,ey,inc,raan,u);

    % Convert from ECI to RTN
    T_ECI2RTN = (rtn2eci_rot(inc,raan,u))'; % rotation matrix from ECI to RTN
    r_analytic_RTN = T_ECI2RTN*r_analytic_ECI;

    % Theorem of Coriolis
    u_dot = sqrt(mu/(sma^3*(1-ecc^2)^3)) * (1+ecc*cos(u))^2;
    w_rtn2eci = [0;0;u_dot]; % rotation of the RTN frame with respect to ECI
    v_analytic_RTN = T_ECI2RTN*v_analytic_ECI - cross(w_rtn2eci,r_analytic_RTN);
    z_analytic(j,:) = [r_analytic_RTN; v_analytic_RTN];

    % Convert numeric position and velocity from ECI to RTN frame
    r_num_ECI = z(j,1:3)';
    v_num_ECI = z(j,4:6)';

    r_num_RTN = (T_ECI2RTN*r_num_ECI);
    v_num_RTN = T_ECI2RTN*v_num_ECI - cross(w_rtn2eci,r_num_RTN);

    % Calculate error between analytical and numerical
    pos_err(j,:) = abs(r_num_RTN-r_analytic_RTN);
    vel_err(j,:) = abs(v_num_RTN-v_analytic_RTN);
end

% Plot
orbit_periods = t / T;

figure(); 
tiledlayout(3,2);

titles = ["Radial (R) Position Error", "Radial (R) Velocity Error", ...
          "Along-Track (T) Position Error", "Along-Track (T) Velocity Error",...
          "Cross-Track (N) Position Error", "Cross-Track (N) Velocity Error"];
ylabels = ["R [m]", "R [m/s]", "T [m]", "T [m/s]", "N [m]", "N [m/s]"];
idxs = [1,1,2,2,3,3];

for i=1:6
    nexttile();
    if rem(i,2) == 0
        plot(orbit_periods, vel_err(:,idxs(i)));
    else
        plot(orbit_periods, pos_err(:,idxs(i)));
    end
    if i == 5 || i == 6
        xlabel("Orbit Periods");
    end
    ylabel(ylabels(i));
    title(titles(i));
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part e
% Osculating orbital elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sma, ecc, inc, raan, aop, ta, u, ecc_bar(3), h_bar(3), h, eps 
osculating = zeros(steps,15);
osculating_j2 = zeros(steps,15);
eps_0 = 1e-10;

for j=1:steps
    % Keplerian orbit
    r = z(j,1:3); v = z(j,4:6);
    % TODO: Update fcns to only give right oe's
    [sma,ecc,inc,raan,aop,ta,u,e_bar,h_bar,eps] = eci2oerad(mu, r, v, eps_0, 1);
    osculating(j,:) = [sma, ecc, rad2deg(inc), rad2deg(raan), rad2deg(aop), rad2deg(ta), rad2deg(u), e_bar, h_bar, norm(h_bar), eps];

    % J2
    r = z_j2(j,1:3); v = z_j2(j,4:6);
    [sma,ecc,inc,raan,aop,ta,u,e_bar,h_bar,eps] = eci2oerad(mu, r, v, eps_0, 1);
    osculating_j2(j,:) = [sma, ecc, rad2deg(inc), rad2deg(raan), rad2deg(aop), rad2deg(ta), rad2deg(u), e_bar, h_bar, norm(h_bar), eps];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part f
% Linear Differential Equations thru OE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical integration
% oe_bar0 = [sma0; ecc0; inc0; raan0; aop0; ta0];
% sma, ecc, ex, ey, inc, raan, ta
eVec0 = (1/ mu ) *(( norm(v0) ^2 - mu / norm(r0)) * r0 - dot ( r0 , v0 )* v0 );
oe_bar0 = [sma0; ecc0; eVec0(1); eVec0(2); inc0; raan0; 0];
z0 = oe_bar0;

% Solve the EOM
% TODO: How to include "n" in du/dt
% is u here mean or true. does it matter?
[t,zoe] = ode45(@(t,z) meanDfq(t,z,rE,mu,J2),tspan,z0,options);

% Eccentricity
zoe(:,8) = zoe(:,3);
zoe(:,9) = zoe(:, 4);
zoe(:,10) = zeros(length(zoe(:, 4)),1);
zoe(:,2) = vecnorm(zoe(:,8:10)');

% Re-arrange so it follows previous set up
zoe(:,3) = rad2deg(zoe(:,5)); % inc
zoe(:,4) = rad2deg(zoe(:,6)); % raan

% Specific angular momentum
zoe(:,11) = NaN; zoe(:,12) = NaN; zoe(:,13) = NaN; 

% zoe(:,6) = f_bar;
% zoe(:,7) = zoe(:,5) + zoe(:,6); % argument of latitude

% Getting angular momentum and energy
p = zoe(:,1) .* (1 - zoe(:,2));
h = sqrt(mu * p);
eps = - mu ./ (2.*zoe(:,1));

zoe(:,14) = h;
zoe(:,15) = eps;
% zoe(:,8) = z(:,6);

% zoe(:,3) = rad2deg(zoe(:,3)); zoe(:,4) = rad2deg(zoe(:,4)); 
% zoe(:,5) = rad2deg(zoe(:,5)); zoe(:,6) = rad2deg(zoe(:,6));
% zoe(:,7) = rad2deg(zoe(:,7));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Parts e, f,g
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
figure();
ylabels = ["a [m]", "e", "i [deg]", "\Omega [deg]",...
           "u [deg]", "||h|| [m^2/s]", "\epsilon [m^2/s^2]"];
cols = [1,2,3,4,7,14,15];

for i=1:7
    subplot(2,4,i); hold on;
    plot(orbit_periods, osculating_j2(:,cols(i)), 'LineWidth',2);
    plot(orbit_periods, osculating(:,cols(i)),'-.', 'LineWidth',2);
    % plot(orbit_periods, zoe(:, cols(i)),'--', 'LineWidth',2)
    xlabel('Orbit Periods');
    ylabel(ylabels(i));
    xlim([0 orb_rev]);
    grid on;
    hold off;
end

subplot(2,4,8);
plot(orbit_periods, NaN);
legend('J2 Perturbed', 'Unperturbed');
% legend('J2 Perturbed', 'Unperturbed', 'Mean classical orbit');
set(gca,'XTick',[])
set(gca,'YTick',[])

% Eccentricity vector
figure();
ylabels = ["e_I", "e_J", "e_K", ... 
           "h_I [m^2/s]", "h_J [m^2/s]", "h_K [m^2/s]"];
cols = [8,9,10,11,12,13];
for i=1:6
    subplot(2,3,i); hold on;
    plot(orbit_periods, osculating_j2(:,cols(i)), 'LineWidth',2);
    plot(orbit_periods, osculating(:,cols(i)), 'LineWidth',2);
    % plot(orbit_periods, zoe(:, cols(i)),'--', 'LineWidth',2)
    xlabel('Orbit Periods');
    ylabel(ylabels(i));
    ytickformat('%.2f')  
    xlim([0 orb_rev]);
    grid on;
    hold off;
end

% Plot u close up
figure(); hold on; 
subplot(3,1,1); grid on;
plot(orbit_periods,zoe(:,7), '--', 'Color', "#EDB120", 'LineWidth',2);
ylabel('u [deg]')
xlabel('Orbit periods')
hold off;

subplot(3,1,2); grid on;
plot(orbit_periods,zoe(:,8), '--', 'Color', "#EDB120", 'LineWidth',2);
ylabel('e_I')
xlabel('Orbit periods')
hold off;

subplot(3,1,3); grid on;
plot(orbit_periods,zoe(:,9), '--', 'Color', "#EDB120", 'LineWidth',2);
ylabel('e_J')
xlabel('Orbit periods')
hold off;













