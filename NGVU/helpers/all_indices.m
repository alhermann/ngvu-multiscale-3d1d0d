function all_indices()
%ALL_INDICES Define index structs 'ind' (states) and 'flu' (pointers).
% Kept as a function to mirror your original code pattern.

    % Make visible to callers without passing explicitly
    % (Most of your code uses these implicitly after calling all_indices).
    clear global ind flu
    global ind flu

    %% ODE state indices (1-based)
    ind.R_k      = 1;
    ind.N_Na_k   = 2;
    ind.N_K_k    = 3;
    ind.N_HCO3_k = 4;
    ind.N_Cl_k   = 5;
    ind.N_Na_s   = 6;
    ind.N_K_s    = 7;
    ind.N_HCO3_s = 8;
    ind.K_p      = 9;
    ind.w_k      = 10;

    ind.Ca_i     = 11;
    ind.s_i      = 12;
    ind.v_i      = 13;
    ind.w_i      = 14;
    ind.I_i      = 15;

    ind.K_i      = 16;

    ind.Ca_j     = 17;
    ind.s_j      = 18;
    ind.v_j      = 19;
    ind.I_j      = 20;

    ind.Mp       = 21;
    ind.AMp      = 22;
    ind.AM       = 23;

    ind.R        = 24;

    % Neuron/NO
    ind.Ca_n       = 25;
    ind.nNOS_act_n = 26;
    ind.NO_n       = 27;
    ind.NO_k       = 28;
    ind.O2_k       = 29;
    ind.NO_i       = 30;
    ind.NO_j       = 31;
    ind.cGMP_i     = 32;
    ind.eNOS_act_j = 33;
    ind.E_b        = 34;
    ind.E_6c       = 35;

    %% Astrocyte pointers (AC)
    flu.R_s       = 1;
    flu.N_Cl_s    = 2;
    flu.Na_k      = 3;
    flu.K_k       = 4;
    flu.HCO3_k    = 5;
    flu.Cl_k      = 6;
    flu.Na_s      = 7;
    flu.K_s       = 8;
    flu.HCO3_s    = 9;
    flu.Cl_s      = 10;

    flu.E_Na_k    = 11;
    flu.E_K_k     = 12;
    flu.E_Cl_k    = 13;
    flu.E_NBC_k   = 14;

    flu.v_k       = 15;

    flu.J_KCC1_k  = 16;
    flu.J_NBC_k   = 17;
    flu.J_NKCC1_k = 18;
    flu.J_NaK_k   = 19;
    flu.J_Na_k    = 20;
    flu.J_K_k     = 21;

    flu.J_BK_k    = 22;
    flu.E_BK_k    = 23;
    flu.w_inf     = 24;
    flu.phi_w     = 25;

    % NO in AC
    flu.p_NO_k    = 26;
    flu.c_NO_k    = 27;
    flu.d_NO_k    = 28;
    flu.tau_ki    = 29;

    %% SMC pointers
    flu.v_coup_i      = 1;
    flu.Ca_coup_i     = 2;
    flu.IP3_coup_i    = 3;
    flu.rho_i         = 4;
    flu.J_IP3_i       = 5;
    flu.J_SRuptake_i  = 6;
    flu.J_CICR_i      = 7;
    flu.J_extrusion_i = 8;
    flu.J_leak_i      = 9;
    flu.J_VOCC_i      = 10;
    flu.J_NaCa_i      = 11;
    flu.J_NaK_i       = 12;
    flu.J_Cl_i        = 13;
    flu.J_K_i         = 14;
    flu.Kactivation_i = 15;
    flu.J_degrad_i    = 16;
    flu.J_stretch_i   = 17;

    flu.v_KIR_i       = 18;
    flu.G_KIR_i       = 19;
    flu.J_KIR_i       = 20;

    flu.M             = 21;
    flu.h_r           = 22;
    flu.E_K_i         = 23;

    % NO / cGMP in SMC
    flu.p_NO_i        = 24;
    flu.c_NO_i        = 25;
    flu.d_NO_i        = 26;
    flu.k_4           = 27;
    flu.E_5c          = 28;
    flu.K_5c          = 29;
    flu.R_cGMP        = 30;
    flu.K_2c          = 31;
    flu.K_act_i       = 32;
    flu.c_w_i         = 33;
    flu.V_max_pde     = 34;
    flu.tau_ij        = 35;

    %% EC pointers
    flu.v_coup_j      = 1;
    flu.Ca_coup_j     = 2;
    flu.IP3_coup_j    = 3;
    flu.rho_j         = 4;
    flu.J_0_j         = 5;
    flu.J_IP3_j       = 6;
    flu.J_ERuptake_j  = 7;
    flu.J_CICR_j      = 8;
    flu.J_extrusion_j = 9;
    flu.J_leak_j      = 10;
    flu.J_cation_j    = 11;
    flu.J_BKCa_j      = 12;
    flu.J_SKCa_j      = 13;
    flu.J_K_j         = 14;
    flu.J_R_j         = 15;
    flu.J_degrad_j    = 16;
    flu.J_stretch_j   = 17;

    % NO in EC
    flu.F_wss         = 18;
    flu.W_wss         = 19;
    flu.tau_wss       = 20;
    flu.p_NO_j        = 21;
    flu.c_NO_j        = 22;
    flu.d_NO_j        = 23;

    %% NE pointers
    flu.w_NR2_A       = 1;
    flu.w_NR2_B       = 2;
    flu.I_Ca          = 3;
    flu.I_Ca_tot      = 4;
    flu.phi_mc        = 5;
    flu.m_c           = 6;
    flu.CaM_n         = 7;   % (not used in all_fluxes; kept for completeness)
    flu.p_NO_n        = 8;
    flu.c_NO_n        = 9;
    flu.d_NO_n        = 10;
    flu.tau_kn        = 11;
end
