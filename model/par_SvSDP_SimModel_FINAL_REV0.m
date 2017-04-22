
%% Parameters for Speed-variable Switched Differential Pump (SvSDP) Cylinder Drive Prototype (with Simplified Check Valves) - Lasse Schmidt, ET/AAU - Date: 04-04-2017
%
clear all
clc
% Simulation Time Step
T_sim = 1/2000;
% Load Pump Coefficients
load('AAAZPF_016_014_011_Pump_Coefficients.mat')
% Load Motion Reference
load('SvSDP_REF_dx125.mat')
% Cylinder Parameters
Dp      = 63E-3;                                % [m] Inner tube diameter (Piston diameter)
Dr      = 35E-3;                                % [m] Outer rod diameter
A_A      = pi/4*Dp^2;                           % [m^2] Piston side Area
A_B      = pi/4*(Dp^2-Dr^2);                    % [m^2] Rod side Area
alpha   = A_B/A_A;                              % [-] Cylinder Area Ratio
VAini   = 0.36e-3;                              % [m^3] Initial Chamber Volume
VBini   = 0.54e-3;                              % [m^3] Initial Chamber Volume
% Mass & Friction Parameters
Ms = 658;                                       % [kg] Mass of Load & Piston
Fcs = 1740.8;                                   % [N] Coulomb Friction
Fss = Fcs+50;                                   % [N] Stribeck Friction
Bs = 6480;                                      % [Ns/m] Viscous Friction Coefficient
gammas =1700;                                   % [-] sign-function smoothing parameters
dx_stribeck = 7e-3;                             % [m/s] Stribeck Friction Characteristic Velocity
% Reservoir (Tank) Pressure
pT = 1e5;                                       % [Pa] Tank Pressure
% Bulk Modulii Parameters
patm = 1.013e5;                                 % [Pa] Athmospheric Pressure
eta_air = 0.7e-2;                               % [%] Percentage Air in Fluid
kappa_air = 1.4;                                % [%] Adiabatic Constant
beta_F = 7500e5;                                % [Pa] Max. Effective Bulk Modulus
% external Load Force
Fext = 20e3;                                    % [N] External Load Force
% Initial Conditions
x_ini   = 50e-3;                                % [m] Initial cylinder position
pB_ini   = 1e5;%30E5;                           % [Pa] Initial B-line pressure
pA_ini   = Fext/A_A+alpha*pB_ini;%30E5;         % [Pa] Initial A-line pressure
% State Limits
x_max   = 700E-3;                               % [m] Maximum cylinder position
x_min   = 0E-3;                                 % [m] Minimum cylinder position
pA_max   = 500E5;                               % [Pa] Maximum A-line pressure
pA_min   = 1E5;                                 % [Pa] Minimum A-line pressure
pB_max   = 500E5;                               % [Pa] Maximum B-line pressure
pB_min   = 1E5;                                 % [Pa] Minimum B-line pressure
% Dynamic Parameters of Servo Drive + PMSM Motor
wn_m = 120*2*pi;                                % [rad/s] Small Signal Bandwidth of Servo Drive
zeta_m = 0.50;                                  % [-] Damping Ratio of Servo Drive
Drive_lim_min = -3000*pi/30;                    % [rad/s] Max. Speed of Servo Drive (Positive Direction)
Drive_lim_max = 3000*pi/30;                     % [rad/s] Max. Speed of Servo Drive (Negative Direction)
Drive_v_lim = 95e3*pi/30;                       % [rad/s] Servo drive slewrate limit
% Dynamic Parameters of 2/2 Valves 
wn_V = 2*pi*21.22;                              % [rad/s] Small Signal Bandwidth of 2/2 Valves
zeta_V = 1;                                     % [-] Damping Ratio of 2/2 Valves
% Parameters for MSK071E-0300-NN Related to Losses
Rs = 0.79;                                      % [ohm] Resistance in PMSM Motor Windings
tau_nomK100 = 28;                               % [Nm] Rated PMSM Motor Shaft Torque
I_nomK100 = 15.8;                               % [A] Current at Rated PMSM Motor Shaft Torque
%% Model Linearization
% Pressures at Linearization Point
pA0 = 25e5;
pB0 = 25e5;
% Piston Position at Linearization Point (Approx. Corresponding to Lowest Eigenfrequency with Constant Inertia Load)
x0 = -(x_max*alpha*A_A+alpha^2*VAini-sqrt(alpha*(x_max*alpha*A_A+alpha*VAini+VBini)^2)+VBini)/((alpha-1)*A_A*alpha);
% Volumes at linearization point
VA_lin = VAini + A_A*x0;                             
VB_lin = VBini + (x_max-x0)*A_B;
% Volume ratio at linearization point
rho = VB_lin/VA_lin;
% Wighting factor H at linearization point
H = VB_lin/(alpha*VA_lin);
% Linear Model Coefficients
Beta_A = 1/(1/beta_F+1/(((1-eta_air)*(patm/(pA0+patm))^(-1/kappa_air)/eta_air+1)*kappa_air*(pA0+patm)));
K_wA = DAp_coeff(1) + DCp_coeff(1);
K_leakA = DAp_coeff(2) + DCp_coeff(2);
Beta_B = 1/(1/beta_F+1/(((1-eta_air)*(patm/(pB0+patm))^(-1/kappa_air)/eta_air+1)*kappa_air*(pB0+patm)));
K_wB = DBp_coeff(1);
K_leakB = DBp_coeff(2);
% Simplified Linear Model Coefficients
Delta_Kw = K_wA - K_wB/alpha;
Lambda_Kw = K_wA + K_wB/H;
K_HpH = alpha*K_leakA+K_leakB/alpha;
K_HpL = H*K_leakA-K_leakB/alpha;
K_LpL = H*K_leakA+K_leakB/H;
K_LpH = alpha*K_leakA-K_leakB/H;
%% Level Pressure Controller
% Level Pressure Transfer Function pH/qH
Gp = tf([Beta_A],[(VAini*((x_max*alpha*A_A+VBini)/(alpha*VAini)+alpha)) Beta_A*K_HpH]);
% Additional Gain of Level Pressure Controller
Kpad = 5;
% Proportional Gain of Level Pressure Controller
Kpp = Kpad*(x_max*alpha*A_A+alpha^2*VAini+VBini)/(alpha*Beta_A);
% Integral Gain of Level Pressure Controller
Kip = Kpad*K_HpH;
% Filter Frequency of Level Pressure Controller 2nd Order Filter
wf = 2*pi*15;
% Damping Ratio of Level Pressure Controller 2nd Order Filter
zetaf = 0.9;%
% Level Pressure Controller
Gppi = tf([Kpp Kip],[1 0])*tf([wf^2],[1 2*zetaf*wf wf^2]);
% Level Pressure Setting
pset = 25e5;
% Margins with Level Pressure Controller
% margin(Gppi*Gp)
%% Motion (Position) Controller
% Gain of Pressure Feedback Loop (Active Damping)
Kad = H/(H+alpha)^2*abs(K_LpL)*10000;
% Gain, Damping Ration % Eigenfrequency of Pressure Feedback Compensated Transformed, Transfer Function xP/qL
wn_ad = sqrt(Beta_A*(A_A^2*(H+alpha)^2+Bs*H*(K_LpL+Kad))/(Ms*H*VA_lin*(H+alpha)));
zeta_ad = (1/2)*H*(Bs*VA_lin*(H+alpha)+Ms*Beta_A*(K_LpL+Kad))*wn_ad/(Beta_A*(A_A^2*(H+alpha)^2+Bs*H*(K_LpL+Kad)));
kn_ad = A_A*(H+alpha)^2/(A_A^2*(H+alpha)^2+Bs*H*(K_LpL+Kad));
% Transfer Function xP/qL (Pressure Feedback Compensated, Transformed Version)
Gx_ad = tf([kn_ad],[1/wn_ad^2 2*zeta_ad/wn_ad 1 0]);
% PI Position Controller Frequency (may be used to adjust e.g. phase margin)
wpi = wn_ad*0.1;
% Desired Gain Margin
eta_gm = 6;
% PI Position Controller Gain with Desired Gain Margin = eta_gm
Kmp = 10^(-eta_gm/20)/(abs((1/2)*kn_ad/(zeta_ad*(2*zeta_ad*wpi-wn_ad)*wpi)));
% PI Position Controller
Gmpi = Kmp*tf([1/wpi 1],[1 0]);
% Margins with PI Position Controller
% margin(Gmpi*Gx_ad)
%% END
