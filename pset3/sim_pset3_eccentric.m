%% We are Close in Eccentric Orbits
% AA 279D Problem Set 3
% Sydney Hsu and Pol Francesch
clear; clc; close all;

initial_conditions_pset3_eccentric;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part a
% Initial TH Orbital Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inertial relative position and velocity
rho_0 = r_init_d - r_init_c; % ECI coordinates
rhodot_0 = v_init_d - v_init_c; % ECI coordinates / time derivative taken in ECI frame

% Relative position and velocity in the RTN frame of the chief
[rho_0_RTN,rhodot_0_RTN] = relative_motion(rho_0,rhodot_0,sma_c,ecc_c,inc_c,raan_c,qns_oe_init_c(6),mu);

state0 = [rho_0_RTN; rhodot_0_RTN];

% Check that this is valid under the Tschauner-Hempel assumptions 
% Small relative separation (rho~0.001r0)
rho_over_r0 = norm(rho_0) / norm(r_init_c);
% Equal semi-major axes
sma_diff = sma_c - oe_init_d1(1);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part b
% YH integration constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ta_c = mean2true(M_c, ecc_c, 1e-8);
stm = dimensional_ya_stm(0, mu, sma_c, ecc_c, ta_c);
K_constants = inv(stm)*state0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Parts c, d
% Yamanaka-Ankersen analytical solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_orbits = 15;
tspan = 0:T/1000:num_orbits*T;
t0 = 0;

% True anomaly as independent variable
state_ya = zeros(length(tspan),6); % in RTN frame

fs = zeros(1,length(tspan));

for j=1:length(tspan)
    % True anomaly
    t = tspan(j);
    M = n*(t - t0);
    f = mean2true(M,ecc_c,1e-8);

    % Relative position and velocity
    state_ya(j,:) = dimensional_ya_stm(t, mu, sma_c, ecc_c, f) * K_constants;
end

% Unbounded motion
orbit_periods = tspan / T;
figure();
subplot(2,1,1);
plot(orbit_periods, state_ya(:,1));
xlabel("Orbital Periods"); ylabel("R (m)"); grid on;

subplot(2,1,2);
plot(orbit_periods, state_ya(:,2));
xlabel("Orbital Periods"); ylabel("T (m)"); grid on;

% Position and velocity RTN plots
figure('Name','YA RTN Position')
plot_rtn_pos(state_ya(:,1:3),'-',1);

figure('Name','YA RTN Velocity')
plot_rtn_vel(state_ya(:,4:6),'-',1); hold on;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part e
% Quasi-nonsingular relative orbital elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From ICs script
roe_init1 = [dsma_1;dlambda_1;dex_1;dey_1;dix_1;diy_1];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part f
% ROE geometric linear mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ex0 = qns_oe_init_c(2); ey0 = qns_oe_init_c(3);

% True anomaly as independent variable
eta = sqrt(1-ecc_c^2);
state_map = zeros(6, length(tspan));

oe = oe_init_c;

for j=1:length(tspan)
    % True anomaly
    t = tspan(j);
    M = n*(t - t0);
    
    if abs(abs(calc_eccentric(M, ecc_c, 1e-5)) - abs(M)) < 1
        E = calc_eccentric(M, ecc_c, 1e-5);
    end

    f = 2 * atan(sqrt((1+ecc_c)/(1-ecc_c)) * tan(E/2));
    % f = mean2true(M, ecc_c, 1e-8);
    % oe(6) = M;
    % qns_oe = singular2qns(oe);
    u = f + aop_c;
    
    % Auxiliary
    k = 1+ex0*cos(u) + ey0*sin(u);
    k_prime = -ex0*sin(u) + ey0*cos(u);
    
    % Coefficients
    bx1 = 1/k + 3/2*k_prime*n/eta^3*t;
    bx2 = -k_prime/eta^3;
    bx3 = 1/eta^3*(ex0*(k-1)/(1+eta)-cos(u));
    bx4 = 1/eta^3*(ey0*(k-1)/(1+eta)-sin(u));
    bx6 = k_prime/eta^3*cot(inc_c);
    by1 = -3/2*k*n/eta^3*t;
    by2 = k/eta^3;
    by3 = 1/eta^2*((1+1/k)*sin(u)+ey0/k+k/eta*(ey0/(1+eta)));
    by4 = -1/eta^2*((1+1/k)*cos(u)+ex0/k+k/eta*(ex0/(1+eta)));
    by6 = (1/k-k/eta^3)*cot(inc_c);
    bz5 = 1/k*sin(u);
    bz6 = -1/k*cos(u);
    bxd1 = k_prime/2+3/2*k^2*(1-k)*n/eta^3*t;
    bxd2 = k^2/eta^3*(k-1);
    bxd3 = k^2/eta^3*(eta*sin(u)+ey0*(k-1)/(1+eta));
    bxd4 = -k^2/eta^3*(eta*cos(u)+ex0*(k-1)/(1+eta));
    bxd6 = -k^2/eta^3*(k-1)*cot(inc_c);
    byd1 = -3/2*k*(1+k*k_prime*n/eta^3*t);
    byd2 = k^2/eta^3*k_prime;
    byd3 = (1+k^2/eta^3)*cos(u)+ex0*k/eta^2*(1+k/eta*(1-k)/(1+eta));
    byd4 = (1+k^2/eta^3)*sin(u)+ey0*k/eta^2*(1+k/eta*(1-k)/(1+eta));
    byd6 = -(1+k^2/eta^3)*k_prime*cot(inc_c);
    bzd5 = cos(u)+ex0;
    bzd6 = sin(u)+ey0;

    % STM
    b_matrix = [bx1, bx2, bx3, bx4, 0, bx6;
                by1, by2, by3, by4, 0, by6;
                0, 0, 0, 0, bz5, bz6;
                bxd1, bxd2, bxd3, bxd4, 0, bxd6;
                byd1, byd2, byd3, byd4, 0, byd6;
                0, 0, 0, 0, bzd5, bzd6];
    state_map(:,j) = [sma_c*eta^2*eye(3), zeros(3,3); zeros(3,3), sma_c*n/eta*eye(3)] * b_matrix * roe_init1;
end

% Plot
state_map = state_map';

% figure();
% rtn_plot_pos([states_c(:,1), state_map(:,1)], ...
%              [states_c(:,2), state_map(:,2)], ...
%              [states_c(:,3), state_map(:,3)]);
% 
% figure();
% rtn_plot_vel([states_c(:,4), state_map(:,4)], ...
%              [states_c(:,5), state_map(:,5)], ...
%              [states_c(:,6), state_map(:,6)]);

figure('Name','RTN position, YA and geometric mapping')
plot_rtn_pos(state_ya(:,1:3),'-',2); hold on;
plot_rtn_pos(state_map(:,1:3),'--',2);
legend('YA','Geometric mapping')

figure('Name','RTN velocity, YA and geometric mapping')
plot_rtn_vel(state_ya(:,4:6),'-',2); hold on;
plot_rtn_vel(state_map(:,4:6),'--',2);
legend('YA','Geometric mapping')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part g
%  Compare integration constants to ROE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_constants_rot = [1 0 0 0 0 0;
       0 cos(aop_c) sin(aop_c) 0 0 0;
       0 -sin(aop_c) cos(aop_c) 0 0 0;
       0 0 0 1 0 0;
       0 0 0 0 cos(aop_c) sin(aop_c);
       0 0 0 0 -sin(aop_c) cos(aop_c)]*K_constants;

icu_2_roe = [1, 0, 0, 0, 0, 0;
             0, -ex0*(eta + 1/(1+eta)), ey0*(eta + 1/(1+eta)), 1, 0, 0;
             0, ex0*ey0, ex0^2-1, -ey0, 0, -ey0*cot(inc_c);
             0, ey0^2-1, ex0*ey0, ex0, 0, ex0*cot(inc_c);
             0, 0, 0, 0, 1, 0;
             0, 0, 0, 0, 0, -1];
roe_from_ic = icu_2_roe * K_constants_rot;

% Difference
roe_ic_diff = roe_from_ic - K_constants

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2 - Part h
%  True relative motion
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
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances

% Solve the EOM
[t,state_num] = ode45(@(t,z) dfq_rel(t,z,mu),tspan,z0,options);

% Plot all together
figure('Name','RTN pos, YA/GM/Numerical');
plot_rtn_pos(state_ya(:,1:3),'-',2); hold on;
plot_rtn_pos(state_map(:,1:3),'--',2); hold on;
plot_rtn_pos(state_num(:,1:3),':',2);

% rtn_plot_pos([state_ya(:,1), state_map(:,1), state_num(:,1)], ...
%              [state_ya(:,2), state_map(:,2), state_num(:,2)], ...
%              [state_ya(:,3), state_map(:,3), state_num(:,3)]);
legend('YA','Geometric mapping','Numerical')

figure('Name','RTN vel, YA/GM/Numerical');
% rtn_plot_vel([state_ya(:,4), state_map(:,4), state_num(:,4)], ...
%              [state_ya(:,5), state_map(:,5), state_num(:,5)], ...
%              [state_ya(:,6), state_map(:,6), state_num(:,6)]);
plot_rtn_vel(state_ya(:,4:6),'-',2); hold on;
plot_rtn_vel(state_map(:,4:6),'--',2); hold on;
plot_rtn_vel(state_num(:,4:6),':',2);
legend('YA','Geometric mapping','Numerical')

% RTN time history combined plot
figure('Name','RTN combined plot')
rtn_plot(orbit_periods, ...
    [state_ya(:,1), state_map(:,1), state_num(:,1)], ...
             [state_ya(:,2), state_map(:,2), state_num(:,2)], ...
             [state_ya(:,3), state_map(:,3), state_num(:,3)], 0, 1);
legend('YA','Geometric mapping','Numerical')

%% Plot differences
% Error from true numerical propagation
ya_error = state_ya - state_num(:,1:6);         % YA solution error
gm_error = state_map - state_num(:,1:6);        % Geometric mapping error

figure();
subplot(3,2,1); hold on; grid on;
plot(orbit_periods, ya_error(:,1)); 
plot(orbit_periods, gm_error(:,1));
ylabel('R error (m)'); hold off;

subplot(3,2,3); hold on; grid on;
plot(orbit_periods, ya_error(:,2)); 
plot(orbit_periods, gm_error(:,2));
ylabel('T error (m)'); hold off;

subplot(3,2,5); hold on; grid on;
plot(orbit_periods, ya_error(:,3)); 
plot(orbit_periods, gm_error(:,3));
ylabel('N error (m)'); xlabel("Orbital Periods"); hold off;

subplot(3,2,2); hold on; grid on;
plot(orbit_periods, ya_error(:,4)); 
plot(orbit_periods, gm_error(:,4));
ylabel('R error (m/s)'); hold off;

subplot(3,2,4); hold on; grid on;
plot(orbit_periods, ya_error(:,5)); 
plot(orbit_periods, gm_error(:,5));
ylabel('T error (m/s)'); hold off;

subplot(3,2,6); hold on; grid on;
plot(orbit_periods, ya_error(:,6)); 
plot(orbit_periods, gm_error(:,6));
legend('YA', 'Geometric mapping');
ylabel('N error (m/s)'); xlabel("Orbital Periods"); hold off;


%% Geometric mapping over time
figure(); 
subplot(2,1,1); hold on;
plot(orbit_periods, state_ya(:,1));
plot(orbit_periods, state_map(:,1),'--');
plot(orbit_periods, state_num(:,1),':');
xlabel("Orbital Periods"); ylabel("R (m)"); grid on;

subplot(2,1,2); hold on;
plot(orbit_periods, state_ya(:,2));
plot(orbit_periods, state_map(:,2),'--');
plot(orbit_periods, state_num(:,2),':');
xlabel("Orbital Periods"); ylabel("T (m)"); grid on;

legend('YA','GA','Num')

%% Numerical solution mapping over time
figure();
subplot(2,1,1);
plot(orbit_periods, state_num(:,1));
xlabel("Orbital Periods"); ylabel("R (m)"); grid on;

subplot(2,1,2);
plot(orbit_periods, state_num(:,2));
xlabel("Orbital Periods"); ylabel("T (m)"); grid on;

