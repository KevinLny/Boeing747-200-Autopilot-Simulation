% Initialization of the Longitudinal Model - Boeing 747-200
clear; clc;

% --- 1. Reference data (Table 4 & 5) ---
S = 510.96;         % Wing surface [m2]
c_bar = 8.32;       % Medium rope [m]
m = 288773.23;      % Mass [kg]
Iyy = 44877574.145; % Pitch inertia [kg.m2]
g = 9.81;           % Gravity [m/s2]
V0 = 205.13;        % VTAS [m/s]
rho = 1.225;        % Air density
q_bar = 13888;      % Dynamic pressure [N/m2]
theta_s = 0.043633; % Stationary pitch angle [rad]
h_m = 6096;         % Altitude [m]
v_tas = 205.13;     

% --- 2. Coefficients Aérodynamiques (Table 6 & 8) ---
CL1 = 0.40; CD1 = 0.0250; CTx1 = 0.0250; 
CDu = 0; CDalpha = 0.20; 
CLu = 0.13; CLalpha = 4.4; CLalpha_dot = 7.0; CLq = 6.6; 
Cmu = 0.013; Cmalpha = -1.0; Cmalpha_dot = -4.0; Cmq = -20.5; 
CL_delta_e = 0.32; Cm_delta_e = -1.30; CD_delta_e = 0; 

% --- 3. Calculation of Dimensional Derivatives (Table 9) ---

Xu = -q_bar*S*(CDu + 2*CD1) / (m*V0);
Xalpha = -q_bar*S*(CDalpha - CL1) / m;
X_delta_e = -q_bar*S*CD_delta_e / m;

Zu = -q_bar*S*(CLu + 2*CL1) / (m*V0);
Zalpha = -q_bar*S*(CLalpha + CD1) / m;
Zalpha_dot = -q_bar*S*c_bar*CLalpha_dot / (2*m*V0);
Zq = -q_bar*S*c_bar*CLq / (2*m*V0);
Z_delta_e = -q_bar*S*CL_delta_e / m;

Mu = q_bar*S*c_bar*(Cmu + 2*Cmu) / (Iyy*V0); 
Malpha = q_bar*S*c_bar*Cmalpha / Iyy;
Malpha_dot = q_bar*S*c_bar^2*Cmalpha_dot / (2*Iyy*V0);
Mq = q_bar*S*c_bar^2*Cmq / (2*Iyy*V0);
M_delta_e = q_bar*S*c_bar*Cm_delta_e / Iyy;

% --- 4. Construction of State Matrices ---
% x = [u; alpha; q; theta]
A = [ Xu,          Xalpha,                0,                -g*cos(theta_s);
      Zu/V0,       Zalpha/V0,             1,                -g*sin(theta_s)/V0;
      Mu+Malpha_dot*Zu/V0, Malpha+Malpha_dot*Zalpha/V0, Mq+Malpha_dot, 0;
      0,           0,                     1,                0 ];

B = [ X_delta_e;
      Z_delta_e/V0;
      M_delta_e + Malpha_dot*Z_delta_e/V0;
      0 ];

C = eye(4);
D = zeros(4,1);

% --- 5. Creation of the System ---
sys_long = ss(A, B, C, D);
sys_tf = tf(sys_long);

% --- Proportional and integral correction of theta ---

Kp_theta = -1.4;
Ki_theta = -0.15;
Kd_theta = -1.1;

fprintf('Modèle longitudinal chargé. Prêt pour Simulink.\n');

% --- Black nichols diagram analysis

s = tf('s');
C_pid = Kp_theta + Ki_theta/s + Kd_theta*s; 
G_avion = ss(A, B, [0,0,0,1], 0);
L = C_pid * G_avion;

%nichols(L);
%grid on;

% --- PID Altitude ---

Kp_h = 1e-3;
Ki_h = 5*0.8e-5;
Kd_h = 4e-2;

% --- Entrée de la commande ---

tau = 5;
cmd_h = 1000;
tau_wind = 5;
