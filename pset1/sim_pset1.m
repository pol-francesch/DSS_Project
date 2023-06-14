%% Orbit Simulation, Review of Astrodynamics
% AA 279D Problem Set 1
% Sydney Hsu and Pol Francesch
clear; clc; close all;

initial_conditions_pset1;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part b
% Initial position and velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[r_init_c,v_init_c] = singular_oe2rv(mu,oe_init_c);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part c
% Numerical integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Numerical integration settings
orb_rev = 5;
stepSize = T/100;
tspan = 0:stepSize:T*orb_rev;
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances
z0 = [r_init_c;v_init_c];

% Solve the EOM
[t,z]       = ode45(@(t,z) dfq(t,z,rE,mu,0),tspan,z0,options);
[t_j2,z_j2] = ode45(@(t,z) dfq(t,z,rE,mu,J2),tspan,z0,options);

set(groot,'defaultAxesFontSize',14);
% Plot
figure('Name', 'Unperturbed'); 
plot_single_orbit_earth(rE,z)

figure('Name', 'J2 Perturbed'); 
plot_single_orbit_earth(rE,z_j2)

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

oe = oe_init_c;

for j = 1:steps
    % Propagate mean anomaly
    M = M0 + n*(t(j)-t0);
    M = wrapTo2Pi(M);
    
    oe(6) = M;
    qns_oe_c = singular2qns(oe);
    u_c = qns_oe_c(6);
    f = mean2true(M, ecc_c, 1e-8);

    % Position & velocity in ECI
    [r_analytic_ECI,v_analytic_ECI] = singular_oe2rv(mu,oe);

    % Convert from ECI to RTN
    T_ECI2RTN = (rtn2eci_rot(inc_c,raan_c,u_c))'; % rotation matrix from ECI to RTN
    r_analytic_RTN = T_ECI2RTN*r_analytic_ECI;

    % Theorem of Coriolis
    f_dot = sqrt(mu/(sma_c^3*(1-ecc_c^2)^3)) * (1+ecc_c*cos(u_c))^2;
    w_rtn2eci = [0;0;f_dot]; % rotation of the RTN frame with respect to ECI
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

%% Plot
orbit_periods = t / T;

figure('Name','RTN error'); 
tiledlayout(3,2);

titles = ["Radial (R) Position Error", "Radial (R) Velocity Error", ...
          "Along-Track (T) Position Error", "Along-Track (T) Velocity Error",...
          "Cross-Track (N) Position Error", "Cross-Track (N) Velocity Error"];
ylabels = ["R (m)", "R (m/s)", "T (m)", "T (m/s)", "N (m)", "N (m/s)"];
idxs = [1,1,2,2,3,3];

for i=1:6
    nexttile();
    if rem(i,2) == 0
        plot(orbit_periods, vel_err(:,idxs(i)));
    else
        plot(orbit_periods, pos_err(:,idxs(i)));
    end
    if i == 5 || i == 6
        xlabel("Orbital Periods");
    end
    ylabel(ylabels(i));
%     xlim([0 10]);
%     title(titles(i));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part e
% Osculating orbital elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sma, ecc, inc, raan, aop, M, u, ecc_bar(3), h_bar(3), h, eps 
osculating = zeros(steps,15);
osculating_j2 = zeros(steps,15);
eps_0 = 1e-10;

for j=1:steps
    % Keplerian orbit
    [oe,eVec,hVec,mechEnergy] = eci2oe(mu, z(j,:));
    qns_oe = singular2qns(oe);
    osculating(j,:) = [oe(1:2)', rad2deg(oe(3:6))', rad2deg(qns_oe(6)), eVec, hVec, norm(hVec), mechEnergy];

    % J2
    [oe,eVec,hVec,mechEnergy] = eci2oe(mu, z_j2(j,:));
    qns_oe = singular2qns(oe);
    osculating_j2(j,:) = [oe(1:2)', rad2deg(oe(3:6))', rad2deg(qns_oe(6)), eVec, hVec, norm(hVec), mechEnergy];
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part f
% Linear Differential Equations thru OE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical integration
% oe_bar0 = [sma0; ecc0; inc0; raan0; aop0; ta0];
% sma, ecc, ex, ey, inc, raan, ta
eVec0 = (1/ mu ) *(( norm(v_init_c) ^2 - mu / norm(r_init_c)) * r_init_c - dot ( r_init_c , v_init_c )* v_init_c );
qns_oe = singular2qns(oe_init_c);
oe_bar0 = [oe_init_c(1:2)'; qns_oe(2:3)'; oe_init_c(3:4)'; qns_oe(6)];
z0 = oe_bar0;

% Solve the EOM
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

% Getting angular momentum and energy
p = zoe(:,1) .* (1 - zoe(:,2));
h = sqrt(mu * p);
eps = - mu ./ (2.*zoe(:,1));

zoe(:,14) = h;
zoe(:,15) = eps;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Parts e, f,g
% Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot
figure('Name','Plots with h and energy');
ylabels = ["$a$ (m)", "$e$", "$i$ (e)", "$\Omega$ (deg)",...
           "$u$ (deg)", "$\big|\big|h\big|\big|$ (m$^2$/s)", "$\epsilon$ (m$^2$/s$^2$)"];
cols = [1,2,3,4,7,14,15];

for i=1:7
    subplot(2,4,i); hold on;
    plot(orbit_periods, osculating_j2(:,cols(i)), 'LineWidth',1);
    plot(orbit_periods, osculating(:,cols(i)),'-.', 'LineWidth',1);
    % plot(orbit_periods, zoe(:, cols(i)),'--', 'LineWidth',1)
    if i > 4
        xlabel('Orbital Periods');
    end
    ylabel(ylabels(i));
    xlim([0 orb_rev]);
    xticks([1:1:orb_rev]);
    grid on;
    hold off;
end

subplot(2,4,8);
plot(orbit_periods, NaN);
lgd = legend('Osculating J2 perturbed', 'Unperturbed');
% lgd = legend('Osculating J2 perturbed', 'Unperturbed', 'Mean J2 perturbed');
fontsize(lgd,10,'points');
set(gca,'XTick',[])
set(gca,'YTick',[])
axis off;
set(gca,'color','none');

% Eccentricity vector
figure('Name','Eccentricity plots');
ylabels = ["$e_I$", "$e_J$", "$e_K$", ... 
           "$h_I$ (m$^2$/s)", "$h_J$ (m$^2$/s)", "$h_K$ (m$^2$/s)"];
cols = [8,9,10,11,12,13];
for i=1:6
    subplot(2,3,i); hold on;
    plot(orbit_periods, osculating_j2(:,cols(i)), 'LineWidth',1);
    plot(orbit_periods, osculating(:,cols(i)), 'LineWidth',1);
    % plot(orbit_periods, zoe(:, cols(i)),'--', 'LineWidth',1)
    if i > 3
        xlabel('Orbital Periods');
    end
    ylabel(ylabels(i));
    ytickformat('%.2f')  
    xlim([0 orb_rev]);
    xticks([1:1:orb_rev]);
    grid on;
    hold off;
end

lgd = legend('Osculating J2 perturbed', 'Unperturbed');
% lgd = legend('Osculating J2 perturbed', 'Unperturbed', 'Mean J2 perturbed');

% Compute variation in hK for unperturbed
hK_std = std(osculating_j2(:,13));

%% Plot u close up
figure('Name','u close up'); hold on; 
subplot(3,1,1); grid on;
plot(orbit_periods,zoe(:,7), '-', 'Color', "#EDB120", 'LineWidth',1);
ylabel('u (deg)')
hold off;

subplot(3,1,2); grid on;
plot(orbit_periods,zoe(:,8), '-', 'Color', "#EDB120", 'LineWidth',1);
ylabel('$e_I$')
hold off;

subplot(3,1,3); grid on;
plot(orbit_periods,zoe(:,9), '-', 'Color', "#EDB120", 'LineWidth',1);
ylabel('$e_J$')
xlabel('Orbital periods')
hold off;



