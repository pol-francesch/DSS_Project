%% Preliminary setup
% Add function paths
addpath('.\pset4\plotting');
addpath('.\plotting');
addpath(genpath('.\src'));

% Constants
J2 = 0.0010826358191967; % Earth's J2 coefficient
mu = 3.986004415e14; % Earth gravitational parameter [m^3/s^2]
rE = 6.378136300e6; % Earth radius [m]

% Plotting
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesFontSize',18);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Integration
options = odeset('RelTol',1e-12,'AbsTol',1e-12); % tolerances

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions - Chief
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial chief absolute singular orbital parameters
sma_c  = 6892.927e3;      % semi-major axis [m] %sma_c = 10000e3;
ecc_c  = 1e-4;            % eccentricity component in I
inc_c  = deg2rad(97.44);  % inclination [rad]
raan_c = deg2rad(270);    % RAAN [rad]
aop_c  = deg2rad(0);      % aop [rad]
M_c    = deg2rad(0);      % mean anomaly [rad]

oe_init_c = [sma_c, ecc_c, inc_c, raan_c, aop_c, M_c]; % combine the above into a single vector
qns_oe_init_c = singular2qns(oe_init_c)'; % chief QNS OE (non-dimensionalized)
oe_init_c_j2 = mean2osc(oe_init_c, 1);

% Position and velocity
[r_init_c,v_init_c] = singular_oe2rv(mu,oe_init_c_j2); % chief position and velocity in inertial frame, J2-perturbed
rv_init_c = [r_init_c;v_init_c];

% Helpful parameters
T = 2*pi*sqrt(sma_c^3/mu); % [sec]
n = sqrt(mu/sma_c^3); % mean motion

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions - ROEs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given ROEs
dsma_1 = 0 / sma_c;
dlambda_1 = 100 / sma_c;
dex_1 = 50 / sma_c;
dey_1 = 100 / sma_c;
dix_1 = 20 / sma_c;
diy_1 = 200 / sma_c;

% Initial conditions - quasi-nonsingular ROEs
roe_init1 = [dsma_1;dlambda_1;dex_1;dey_1;dix_1;diy_1]; % non-dimensionalized

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions - Deputy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial pos/vel of deputy
oe_init_d1 = roe2singular(oe_init_c,roe_init1);
qns_oe_init_d1 = singular2qns(oe_init_d1);

oe_init_d1_j2 = mean2osc(oe_init_d1, 1);

% Position and velocity
[r_init_d,v_init_d] = singular_oe2rv(mu,oe_init_d1_j2); % deputy position and velocity in inertial frame, J2 perturbed
rv_init_d = [r_init_d;v_init_d];
