% Given Transfer Function G
numG = [0,0,292.423680000000];
denG = [1,14.1176111111111,1.25920326222223];
G = tf(numG, denG);

% Define the Laplace variable 's'
s = zpk('s');

% Define new W1 and W3
wS = createSweight(1,1e-4,1.5,1);
wT = createTweight(10,1e-4,1.5,1);

% Augment the plant with weights
P = augw(G, wS, [], wT);

% Perform H-infinity synthesis
gmin = 0.1;  % Minimum gamma value
gmax = 10;   % Maximum gamma value
tol = 0.01;  % Tolerance for gamma convergence

% Perform H-infinity synthesis
[K, CL, gamma] = hinfsyn(P, 1, 1, [gmin gmax]);

% Display the final gamma value
gamma


% Convert the controller K to transfer function form and minimize it
sys = ss(K.A, K.B, K.C, K.D);
sys = minreal(tf(sys));

% Display the resulting controller transfer function
sys

%% Results analyzed

looptransfer = loopsens(G,K);
L = looptransfer.Lo;
T = looptransfer.To;
I = eye(size(L));

figure(1)
omega = logspace(-1,3,100);
sigma(I+L,'b-',wS/gamma,'r--',T,'b-.',gamma/wT,'r.',omega)
grid
legend('1/\sigma(S) performance', ...
'\sigma(wS) performance bound', ...
'\sigma(T) robustness', ...
'\sigma(1/wT) robustness bound')
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times')
figure(2)
omega = logspace(-1,3,100);
sigma(L,'b-',wS/gamma,'r--',gamma/wT,'r.',omega)
grid
legend('\sigma(L)','\sigma(wS) performance bound', ...
'\sigma(1/wT) robustness bound')
set(gca,'FontSize',9,'Fontsize',14,'FontName','Times')

%% Hinf COntrol with Models
% numG = 6.548904113729382e+02;
% denG = [1,2.042422007272920e+02,2.882160103432319e+03,2.571536346057749e+02];
% G = tf(numG,denG);
numG = [0,0,292.423];
denG = [1,14.117,1.259];
G = tf(numG, denG);
nuM = 1;
deM = [2.5 1];
M = tf(nuM,deM);
tol = 1e-2;
nuWp = [1 0.1];
dnWp =[0.01 10];
Wp = tf(nuWp,dnWp);
nuWu = [0.01 0.01];
dnWu = [1 1];
Wu = tf(nuWu,dnWu);

%% Plant augmentation
% Assign input and output names for each system
systemnames = 'G M Wp Wu';
inputvar = '[ ref; dist; control ]';
outputvar = '[ Wp; Wu; ref-G-dist ]';
input_to_G = '[ control ]';
input_to_M = '[ ref ]';
input_to_Wp = '[ G+dist-M ]';
input_to_Wu = '[ control ]';
sys_ic = sysic;

%% Controller Design
nmeas = 1;
ncont = 1;
gmin = 0.1;
gmax = 10;
tol = 0.001;
[K, clp, gamma] = hinfsyn(sys_ic, nmeas, ncont, [gmin, gmax]);
get(K)


%% Hinf design 
% generate_hinf.m
lp_ic = lft(sys_ic,K);
omega = logspace(-2,2,100);
clp_g = ufrd(clp,omega);
opt = robopt('Display','on');
[stabmarg,destabu,report,info] = robuststab(clp_g,opt);
report
semilogx(info.MussvBnds(1,1),info.MussvBnds(1,2))


%% Create state space model from Simulink
sys_G=ss(G); %% Get stace-space model from plantM simulink
[Am,Bm,Cm,Dm] = ssdata(sys_G);
M_simulink=ltisys(Am,Bm,Cm,Dm); %% convert state-space to SYS model M

%% add two unstable poles previously removed from Wt.
% always need 2 pole removales because Wt has 2 poles
M_filtered=sderiv(M_simulink,2,[1/abs(Wt.zpk.value.z{1}(1)) 1 ]);
M_filtered=sderiv(M_filtered,2,[1/abs(Wt.zpk.value.z{1}(2)) 1 ]);

% M_min = ss(M_filtered); % minreal of SYS model after adding poles

%% Find set of solutions to matrix inequalities
% [gopt,Cmod]=hinflmi(M_filtered,[1 1]); % simple
[gopt,Cmod]=hinflmi(M_filtered,[1,1],0,1e-6,[0 0 0]); % with more args

[Ac,Bc,Cc,Dc]=ltiss(Cmod);
[Gc.tf.num,Gc.tf.den] = ss2tf(Ac,Bc,Cc,Dc);
Gc.tf.value=tf(Gc.tf.num,Gc.tf.den);
Gc.zpk.value=zpk(Gc.tf.value);

%% Modify controller if necessary
Gc.gainmod=1;

% Change low frequency poles to s^mu
Gc.mod.value=Gc.tf.value*Gc.gainmod*(s+0.01)^sys.mu/s^sys.mu; % does have poles at origin
% Gcmod = Gc_zpk*(s-Gc_zpk.p{1}(2))*(s-Gc_zpk.p{1}(3))/(s+abs(Gc_zpk.p{1}(2))) / (s+abs(Gc_zpk.p{1}(3)))
% Gcmod=minreal(Gcmod,1e-4)

% Gcmod=Gcmod*(s+9.859)^sys_mu/s^sys_mu % does have poles at origin

% Gcmod=Gc*Gc_gainmod; % doesnt have poles at origin
Gc.mod.value=minreal(Gc.mod.value,1e-4);
Gc.mod.zpk=zpk(Gc.mod.value);

% Gcmod=(0.4353*s^3 + 240.0*s^2 + 11222.0*s + 687.9)/(s^3 + 7218.0*s^2 + 141990.0*s + 698550.0)


%% Generate complementary functions

% Generate nominal loop and sensitivity functions
L.nominal.value=Gc.mod.value*Ga*Gp.nominal.tf*Gs*Gf;
T.nominal.value=L.nominal.value/(1+L.nominal.value);
S.nominal.value=1-T.nominal.value;

% Calculate controller DC gain
Kc=dcgain((s^sys.mu)*Gc.mod.value);

S.nominal.star=S.nominal.value/s^(sys.h); % Remove poles and zeros at s=0
S.nominal.star=minreal(S.nominal.star);% Cancel poles and zeros at s=0
