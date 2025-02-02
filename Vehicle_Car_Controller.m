%% Initialization
clear
load_system('Plant_Model');
freq_range = logspace(-2, 4, 501); % useful for bode & sigma plots.

%% Linearize model
% Linearize your plant model 'Plant_Model'
[a_mat, b_mat, c_mat, d_mat] = linmod('Plant_Model');
P = ss(a_mat, b_mat, c_mat, d_mat);
P.OutputName = {'VehicleSpeed'};  % Replace with actual output names
clear a_mat b_mat c_mat d_mat

figure(1)
clf
sigma(P, 'b-', freq_range)
ylim([-70 70]);
grid on


%% Design the closed-loop response
% Without imaginary poles, directly design the closed-loop target.
s = zpk('s');

% Assuming simple dynamics without imaginary poles:
% Use a simple low-pass filter for closed-loop design.
om_filt = 8;  % Experimentally chosen
M = zpk([], [-om_filt -om_filt -om_filt -om_filt], 1);  
M = M / (dcgain(M) * dcgain(P(1,:)));

R_ideal = P * M;
R_ideal = minreal(R_ideal, [], false);

figure(1)
clf
subplot(121)
step(P, 'b-', 100);
title('Plant Step Response')
subplot(122)
step(R_ideal, 'm-', 100)
title('Target Step Response')

%% Design weights for H-infinity loop-shaping
om_c = bandwidth(P(1,:) * M);  % Same as closed-loop target, can be higher

% Initial attempt at weights - get the overall shape
Wo = 1;
Wi = 1/s;
Pw = Wo * P * Wi;
Pw = minreal(Pw, [], false);

% Adjust Wi to set the closed-loop bandwidth as desired
sg = max(sigma(Pw, om_c)); 
Wi = Wi / sg;
Pw = Pw / sg;

figure(1)
clf
sigma(P, 'b-', Pw, 'g-', freq_range)
ylim([-70 70]);
grid on
legend('P', 'Pw', 'Location', 'NorthEast')

%% Synthesize a controller using ncfsyn
[Cinf, ~, gam] = ncfsyn(Pw);
Cinf = -Cinf;
perf_margin = 1 / gam 

% Form the overall controller
C = Wi * Cinf * Wo;

% Work out the overall loop gain
L = P * C;

figure(1)
clf
sigma(P, 'b', Pw, 'g', L, 'm', freq_range)
ylim([-70 70]);
grid on


