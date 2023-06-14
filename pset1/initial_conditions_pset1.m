%% Preliminary setup
% Add function paths
addpath('.\pset1\plotting');
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

% Helpful parameters
T = 2*pi*sqrt(sma_c^3/mu); % [sec]
n = sqrt(mu/sma_c^3); % mean motion