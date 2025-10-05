%% ===================== Global / Physical constants ===================== %%
Farad    = 9.65e4;     % [C mol^-1]  Faraday constant
R_gas    = 8.315;      % [J mol^-1 K^-1]
Temp     = 300;        % [K]
unitcon  = 1e3;        % [-]         Convenience unit conversion factor

%% ===================== Ostby astrocyte / ECS constants ================= %%
L_p       = 2.1e-9;    % [m uM^-1 s^-1]
R_tot     = 8.79e-8;   % [m]
X_k       = 12.41e-3;  % [uM·m]
z_Na      = 1;         % [-]
z_K       = 1;         % [-]
z_Cl      = -1;        % [-]
z_NBC     = -1;        % [-]
g_K_k     = 40;        % [ohm^-1 m^-2]
g_KCC1_k  = 1e-2;      % [ohm^-1 m^-2]
g_NBC_k   = 7.57e-1;   % [ohm^-1 m^-2]
g_Cl_k    = 8.797e-1;  % [ohm^-1 m^-2]
g_NKCC1_k = 5.54e-2;   % [ohm^-1 m^-2]
g_Na_k    = 1.314;     % [ohm^-1 m^-2]
J_NaK_max = 1.42e-3;   % [uM·m s^-1]
K_Na_k    = 1.0e4;     % [uM]
K_K_s     = 1.5e3;     % [uM]
k_C       = 7.35e-5;   % [uM s^-1]

%% =========================== BK channel (AC) =========================== %%
A_ef_k   = 3.7e-9;     % [m^2]  Endfoot area
v_6      = 22e-3;      % [V]
v_4      = 14.5e-3;    % [V]
psi_w    = 2.664;      % [s^-1]
G_BK_k   = 4.3e3;      % [pS]
g_BK_k   = G_BK_k * 1e-12 / A_ef_k;  % [ohm^-1 m^-2]
VR_pa    = 1e-3;       % [-]     Vol. ratio PVS/AC
VR_ps    = 1e-3;       % [-]     Vol. ratio PVS/SMC

%% ===================== Smooth muscle cell (SMC) ======================== %%
% Filosa-fit parameters
F_il = 7.5e2;    % [-]
z_1  = 4.5;      % [-]
z_2  = -1.12e2;  % [-]
z_3  = 4.2e-1;   % [-]
z_4  = -1.26e1;  % [-]
z_5  = -7.4e-2;  % [-]

% Koenigsberger et al. (units in comments)
Fmax_i = 0.23;        % [uM/s]
Kr_i   = 1;           % [uM]
G_Ca   = 0.00129;     % [uM/mV/s]
v_Ca1  = 100;         % [mV]
v_Ca2  = -35;         % [mV]
R_Ca   = 8.5;         % [mV]
G_NaCa = 0.00316;     % [uM/mV/s]
c_NaCa = 0.5;         % [uM]
v_NaCa = -30;         % [mV]
B_i    = 2.025;       % [-]
cb_i   = 1;           % [uM]
C_i    = 55;          % [-]
sc_i   = 2;           % [uM]
cc_i   = 0.9;         % [uM]
D_i    = 0.24;        % [-]
vd_i   = -100;        % [mV]
Rd_i   = 250;         % [mV]
L_i    = 0.025;       % [1/s]
gam    = 1970;        % [mV uM^-1]  Membrane potential scaling
F_NaK  = 0.0432;      % [-]
G_Cl   = 0.00134;     % [1/mV/s]
v_Cl   = -25;         % [mV]
G_K    = 0.00446;     % [1/mV/s]
vK_i   = -94;         % [mV]
lab    = 45;          % [1/s]
c_w    = 0;           % [-]
bet    = 0.13;        % [-]
v_Ca3  = -27;         % [mV]
R_K    = 12;          % [mV]
k_i    = 0.1;         % [1/s]
K_d    = 1;           % [uM]
B_T    = 100;         % [uM]
Kinf_i = 1e5;         % [uM]  (100 mM)

% Stretch / mechanoelectric coupling
G_stretch = 0.0061;   % [uM mV^-1 s^-1]
P_str     = 30;       % [-]
Esac      = -18;      % [mV]
alpha1    = 0.0074;   % [-]
sig0      = 500;      % [-]

%% ======================= Endothelial cell (EC) ========================= %%
Fmax_j = 0.23;   % [uM/s]
Kr_j   = 1;      % [uM]
B_j    = 0.5;    % [-]
cb_j   = 1;      % [uM]
C_j    = 5;      % [-]
sc_j   = 2;      % [uM]
cc_j   = 0.9;    % [uM]
D_j    = 0.24;   % [-]
L_j    = 0.025;  % [1/s]

G_cat  = 0.66e-3;  % [-]
E_Ca   = 50;       % [mV]
m3cat  = -0.18;    % [-]
m4cat  = 0.37;     % [-]
J0_j   = 0.029;    % [uM/s]  Constant Ca influx (EC)

C_m    = 25.8;     % [-]
G_tot  = 6927;     % [-]
vK_j   = -80;      % [mV]

a1 = 53.3; a2 = 53.3; b = -80.8; c = -6.4;
m3b = 1.32e-3; m4b = 0.3;
m3s = -0.28;   m4s = 0.389;

G_R    = 955;      % [-]
v_rest = -31.1;    % [mV]
k_j    = 0.1;      % [1/s]

%% ========== Gap coupling / IP3 coupling (by scenario CASE) ============= %%
global CASE
switch CASE
    case 1
        g_hat = 50;   p_hat = 0;    p_hatIP3 = 0.05;
    case 2
        g_hat = 5;    p_hat = 0.05; p_hatIP3 = 0.05;
    case 3
        g_hat = 0;    p_hat = 0;    p_hatIP3 = 0.05;
    case 0
        g_hat = 0;    p_hat = 0;    p_hatIP3 = 0;
    otherwise
        % Default to CASE 2 if unset
        g_hat = 5;    p_hat = 0.05; p_hatIP3 = 0.05;
end

%% ===================== Myosin cross-bridge model ======================== %%
K2_c      = 0.5;   % [1/s]
K3_c      = 0.4;   % [1/s]
K4_c      = 0.1;   % [1/s]
K5_c      = 0.5;   % [1/s]
K7_c      = 0.1;   % [1/s]
gam_cross = 17;    % [-]

%% ============================ Kelvin–Voigt ============================== %%
P_r     = 4000;     % [-]
rb_r    = 20e-6;    % [m]
h0_r    = 3e-6;     % [m]
R0pas_r = 2e-6;     % [m]
R0act_r = 12e-6;    % [m]
Etot_r  = 233e3;    % [Pa]
Eact_r  = 23.3e3;   % [Pa]
Epas_r  = 6.6e3;    % [Pa]
nu_r    = 1e4;      % [-]

%% ============================ NO Simulation ============================= %%
% Diffusion / geometry
D_c_NO = 3300;   % [um^2/s]

% Astrocyte / NE interface / neuronal Ca dynamics
Glu_max    = 1846;     % [uM]
V_spine    = 8e-8;     % [nL]
kappa_ex   = 1.6e3;    % [1/s]
Ca_rest    = 0.1;      % [uM]
lam_buf    = 20;       % [-]
V_max_nNOS = 25e-3;    % [uM/s]
K_m_nNOS   = 9.27e-2;  % [-]
mu_deact_n = 0.0167;   % [1/s]

K_m_A   = 650;     % [uM]
K_m_B   = 2800;    % [uM]
v_n     = -0.04;   % [V]
G_M     = 4.6e4;   % [fS]
P_Ca_P_M = 3.6;    % [-]
Ca_ex   = 2e3;     % [uM]
M       = 1.3e5;   % [uM]
alpha_v = -80;     % [V^-1]
beta_v  = 0.02;    % [V]
n_NR2_A = 0.63;    % [-]
n_NR2_B = 11;      % [-]

Q_1 = 1.9e5; Q_2 = 2.1e5; Q_3 = 0.4e5; Q_4 = 0.26e5;  % [uM^-1]
V_max_NO_n = 4.22;      % [1/s]
O2_n       = 200;       % [uM]
K_m_O2_n   = 243;       % [uM]
LArg_n     = 100;       % [uM]
K_m_LArg_n = 1.5;       % [uM]
k_O2_n     = 9.6e-6;    % [uM^-2 s^-1]
x_nk       = 25;        % [um]

% Astrocyte / SMC interface
k_O2_k = 9.6e-6;  % [uM^-2 s^-1]
x_ki   = 25;      % [um]

% SMC NO/cGMP pathway
k_m1      = 100;       % [1/s]
k_1       = 2e3;       % [uM^-1 s^-1]
k_2       = 0.1;       % [1/s]
k_3       = 3;         % [uM^-1 s^-1]
V_max_sGC = 0.852;     % [uM/s]
K_m_pde   = 2;         % [uM]
k_dno     = 0.01;      % [1/s]
C_4       = 0.011;     % [s^-1 uM^-2]
m_4       = 2;         % [-]
K_m_mlpc  = 5.5;       % [uM]
delta_i   = 58.1395;   % [-]
k_mlpc_b  = 0.0086;    % [1/s]
k_mlpc_c  = 0.0327;    % [1/s]
alpha_i   = 1e14;      % [uM^-1]
beta_i    = 0.13;      % [uM^2]
gam_i     = -3;        % [uM^-1]
epsilon_i = 1;         % [uM^-1]
v_Ca3_i   = -27;       % [mV]
R_K_i     = 12;        % [mV]
k_pde     = 0.0195;    % [1/s]
x_ij      = 3.75;      % [um]

% Endothelium NO production / WSS pathway
gam_eNOS  = 0.1;     % [-]
mu_deact_j = 0.0167; % [1/s]
K_dis     = 0.09;    % [uM s^-1]
K_m_eNOS  = 0.45;    % [uM]
g_max     = 0.06;    % [uM s^-1]
alpha_wss = 2;       % [-]
W_0       = 1.4;     % [Pa^-1]
delta_wss = 2.86;    % [Pa]
V_max_NO_j = 1.22;   % [1/s]
O2_j      = 200;     % [uM]
K_m_O2_j  = 7.7;     % [uM]
LArg_j    = 100;     % [uM]
K_m_LArg_j = 1.5;    % [uM]
DeltaPL   = 9.1e4;   % [Pa m^-1]
k_O2_j    = 9.6e-6;  % [uM^-2 s^-1]
