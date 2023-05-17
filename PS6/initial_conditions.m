%% Reconfiguration simulation from mode 1 to mode 2 using least squares
% AA 279D Problem Set 5
% Sydney Hsu and Pol Francesch

%% Preliminary setup
% Add function paths
addpath('.\plots');
addpath('.\mean_osc');
addpath('.\mappings');
addpath('.\dfqs');
addpath('.\frames');

% Constants
J2 = 0.0010826358191967; % Earth's J2 coefficient
mu = 3.986004415e14; % Earth gravitational parameter [m^3/s^2]
rE = 6.378136300e6; % Earth radius [m]

% Plotting
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesFontSize',18);

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
% Mode 1: Preliminary DEM generation offsets (phase C1)
dsma_1 = 0;
dlambda_1 = 0;
de_1 = 260 / sma_c;
di_1 = 222 / sma_c;
phase_diff_1 = 200; % designed phase difference % TN plane not relevant as it is safe in the other planes; along-track is the greatest uncertainty
theta_1 = deg2rad(90); % phase angle of inclination

% Calculate deputy initial orbit parameters
phi_1 = theta_1 - deg2rad(phase_diff_1); % phase angle of eccentricity
dex_1 = de_1*cos(phi_1);  dey_1 = de_1*sin(phi_1);
dix_1 = di_1*round(cos(theta_1),12); diy_1 = di_1*sin(theta_1);
% Initial conditions - quasi-nonsingular ROEs
roe_init1 = [dsma_1;dlambda_1;dex_1;dey_1;dix_1;diy_1]; % non-dimensionalized

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions - Deputy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial pos/vel of deputy
% qns_oe_init_d1 = roe2qns(qns_oe_init_c,roe_init1); % QNS orbital elements of the deputy
% oe_init_d1 = qns2singular(qns_oe_init_d1); % singular orbital elements of the deputy [sma, ecc, inc, raan, aop, ta]
oe_init_d1 = roe2singular(oe_init_c,roe_init1);
qns_oe_init_d1 = singular2qns(oe_init_d1);

oe_init_d1_j2 = mean2osc(oe_init_d1, 1);

% Position and velocity
[r_init_d,v_init_d] = singular_oe2rv(mu,oe_init_d1_j2); % deputy position and velocity in inertial frame, J2 perturbed
rv_init_d = [r_init_d;v_init_d];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify ICs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% verify_oe_c = rv2singular_oe(mu, r_init_c, v_init_c);
% verify_oe_d = rv2singular_oe(mu, r_init_d, v_init_d);

true_roe = roe_init1'*sma_c;

% Given a set of Singular OE for the chief and QNS ROE:
% Convert to QNS OE of chief
% Convert to QNS OE of deputy
% Convert back to ROEs
% Sucess!
verify_roe0 = oes2roe(qns_oe_init_c, qns_oe_init_d1)*sma_c;

% Given a set of singular OE of the chief and QNS ROE:
% Convert to QNS OE of chief
% Convert to QNS OE of deputy
% Convert to singular OE of deputy
% Convert back to ROE
% Broken (either qns -> sing OR sing oe -> qns roe)
verify_roe1 = singular_oe2roe(oe_init_c,oe_init_d1)*sma_c;

% Given a set of singular OE of the chief and QNS ROE:
% Convert to QNS OE of chief
% Convert to QNS OE of deputy
% Convert to singular OE of deputy
% Convert from mean to osculating singular elements
% Convert to RV in ECI
% Convert back to singular OE osc
% Convert to QNS OE osc
% Convert to ROE osc
% Broken
[verify_oe_c,~,~,~] = eci2oe(mu, rv_init_c);
[verify_oe_d,~,~,~] = eci2oe(mu, rv_init_d);
verify_qns_c = singular2qns(verify_oe_c);
verify_qns_d = singular2qns(verify_oe_d);

verify_roe2 = oes2roe(verify_qns_c, verify_qns_d)*sma_c;

% Given a set of singular OE of the chief and QNS ROE:
% Convert to QNS OE of chief
% Convert to QNS OE of deputy
% Convert to singular OE of deputy
% Convert from mean to osculating singular elements
% Convert to QNS OE
% Convert to QNS ROE
% dlambda broken. Everything else is reasonable for osculating ROE
% verify_qns_c2 = singular2qns(oe_init_c_j2);
% verify_qns_d2 = singular2qns(oe_init_d1_j2);
% 
% verify_roe3 = oes2roe(verify_qns_c2, verify_qns_d2)*sma_c;
% 
% % 
% verify_roe4 = singular_oe2roe(oe_init_c_j2,oe_init_d1_j2)*sma_c;





