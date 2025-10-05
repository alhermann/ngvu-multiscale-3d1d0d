function dy = NVC(time, state, ca, G_syn, tmax, t_start, t_end, varargin)
% NVC RHS of the merged neurovascular coupling ODE system.
%
% Inputs:
%   time   : current time [s]
%   state  : state vector (see all_indices)
%   ca     : astrocyte Ca trace (vector or scalar); used by crossbridge kinetics
%   G_syn  : synaptic glutamate trace (vector or scalar); drives NE block
%   tmax   : total sim time [s] (used to convert trace length to dt)
%   t_start, t_end : stimulus window (passed to all_fluxes)
%
% Output:
%   dy     : time-derivative vector

    all_constants();
    all_indices();

    dy = zeros(size(state));

    %% ----- Synaptic drive indexing (Glu) -----
    if numel(G_syn) == 1
        G_syn_index = 1;
    else
        G_syn_data_dt = (tmax * 1000) / length(G_syn);  % ms
        G_syn_index = round(time * 1000 / G_syn_data_dt);
        G_syn_index = max(1, min(G_syn_index, length(G_syn)));
    end
    Glu_sc = G_syn(G_syn_index) * 1e3;  % mM → uM

    %% ----- Fluxes / algebraic relations -----
    [NE, AC, SMC, EC] = all_fluxes(time, state, t_start, t_end, Glu_sc);

    %% ----- Neuron NO-pathway -----
    CaM_n              = (state(ind.Ca_n)) / NE(flu.m_c);   % [uM] CaM-Ca
    dy(ind.Ca_n)       = NE(flu.I_Ca_tot) / (2 * Farad * V_spine) ...
                       - kappa_ex * (state(ind.Ca_n) - Ca_rest) / (1 + lam_buf);
    dy(ind.nNOS_act_n) = (V_max_nNOS * CaM_n) / (K_m_nNOS + CaM_n) - mu_deact_n * state(ind.nNOS_act_n);
    dy(ind.NO_n)       = NE(flu.p_NO_n) - NE(flu.c_NO_n) + NE(flu.d_NO_n);

    %% ----- Astrocyte -----
    dy(ind.R_k)      = L_p * ( AC(flu.Na_k) + AC(flu.K_k) + AC(flu.Cl_k) + AC(flu.HCO3_k) ...
                             - AC(flu.Na_s) - AC(flu.K_s) - AC(flu.Cl_s) - AC(flu.HCO3_s) ...
                             + X_k / state(ind.R_k) );

    dy(ind.N_Na_k)   = - AC(flu.J_Na_k)  - 3*AC(flu.J_NaK_k) + AC(flu.J_NKCC1_k) + AC(flu.J_NBC_k);
    dy(ind.N_K_k)    = - AC(flu.J_K_k)   + 2*AC(flu.J_NaK_k) + AC(flu.J_NKCC1_k) + AC(flu.J_KCC1_k) ...
                       - AC(flu.J_BK_k);
    dy(ind.N_HCO3_k) = 2 * AC(flu.J_NBC_k);
    dy(ind.N_Cl_k)   = dy(ind.N_Na_k) + dy(ind.N_K_k) - dy(ind.N_HCO3_k);

    dy(ind.N_Na_s)   = - 0.0 - dy(ind.N_Na_k);
    dy(ind.N_K_s)    =   0.0 - dy(ind.N_K_k);
    dy(ind.N_HCO3_s) = - dy(ind.N_HCO3_k);

    dy(ind.K_p)      = AC(flu.J_BK_k) / (VR_pa * state(ind.R_k)) + SMC(flu.J_KIR_i) / (VR_ps);
    dy(ind.w_k)      = AC(flu.phi_w) * (AC(flu.w_inf) - state(ind.w_k));

    % AC NO
    dy(ind.NO_k)     = AC(flu.p_NO_k) - AC(flu.c_NO_k) + AC(flu.d_NO_k);

    %% ----- SMC -----
    dy(ind.Ca_i) = SMC(flu.Ca_coup_i) + SMC(flu.rho_i) * ( SMC(flu.J_CICR_i) + SMC(flu.J_IP3_i) ...
                 + SMC(flu.J_leak_i) - SMC(flu.J_SRuptake_i) - SMC(flu.J_extrusion_i) ...
                 - SMC(flu.J_VOCC_i) + SMC(flu.J_NaCa_i) + 0.1*SMC(flu.J_stretch_i) );

    dy(ind.s_i)  = - SMC(flu.J_CICR_i) - SMC(flu.J_leak_i) + SMC(flu.J_SRuptake_i);
    dy(ind.v_i)  = SMC(flu.v_coup_i) + gam * ( - SMC(flu.J_NaK_i) - SMC(flu.J_Cl_i) ...
                   - 2*SMC(flu.J_VOCC_i) - SMC(flu.J_NaCa_i) - SMC(flu.J_K_i) ...
                   - SMC(flu.J_stretch_i) - SMC(flu.J_KIR_i) );
    dy(ind.w_i)  = lab * (SMC(flu.Kactivation_i) - state(ind.w_i));
    dy(ind.I_i)  = SMC(flu.IP3_coup_i) - SMC(flu.J_degrad_i);

    dy(ind.K_i)  = - SMC(flu.J_KIR_i) - SMC(flu.J_K_i) + SMC(flu.J_NaK_i);

    % NO / cGMP
    dy(ind.NO_i)     = SMC(flu.p_NO_i) - SMC(flu.c_NO_i) + SMC(flu.d_NO_i);
    dy(ind.E_b)      = -k_1 * state(ind.E_b) * state(ind.NO_i) + k_m1 * state(ind.E_6c) + SMC(flu.k_4) * SMC(flu.E_5c);
    dy(ind.E_6c)     =  k_1 * state(ind.E_b) * state(ind.NO_i) - (k_m1 + k_2) * state(ind.E_6c) - k_3 * state(ind.E_6c) * state(ind.NO_i);
    dy(ind.cGMP_i)   =  V_max_sGC * SMC(flu.E_5c) - (SMC(flu.V_max_pde) * state(ind.cGMP_i)) / (K_m_pde * state(ind.cGMP_i));

    %% ----- EC -----
    dy(ind.Ca_j) = EC(flu.Ca_coup_j) + EC(flu.rho_j) * ( EC(flu.J_IP3_j) - EC(flu.J_ERuptake_j) ...
                 + EC(flu.J_CICR_j) - EC(flu.J_extrusion_j) + EC(flu.J_leak_j) ...
                 + EC(flu.J_cation_j) + EC(flu.J_0_j) + EC(flu.J_stretch_j) );
    dy(ind.s_j)  = EC(flu.J_ERuptake_j) - EC(flu.J_CICR_j) - EC(flu.J_leak_j);
    dy(ind.v_j)  = - ( EC(flu.J_K_j) + EC(flu.J_R_j) ) / C_m + EC(flu.v_coup_j);
    dy(ind.I_j)  = EC(flu.IP3_coup_j) + J_PLC - EC(flu.J_degrad_j);

    % EC NO
    dy(ind.eNOS_act_j) = gam_eNOS * (K_dis * state(ind.Ca_j)) / (K_m_eNOS + state(ind.Ca_j)) ...
                       + (1 - gam_eNOS) * g_max * EC(flu.F_wss) - mu_deact_j * state(ind.eNOS_act_j);
    dy(ind.NO_j)       = EC(flu.p_NO_j) - EC(flu.c_NO_j) + EC(flu.d_NO_j);

    %% ----- Crossbridge (depends on astrocyte Ca) -----
    if numel(ca) == 1
        ca_index = 1;
    else
        ca_data_dt = (tmax * 1000) / length(ca);  % ms
        ca_index   = round(time * 1000 / ca_data_dt);
        ca_index   = max(1, min(ca_index, length(ca)));
    end

    K1_c = gam_cross * (ca(ca_index)*1e-3)^3;
    K6_c = K1_c;

    % NO-mediated dilation
    K2_c = delta_i * (k_mlpc_b + k_mlpc_c * SMC(flu.R_cGMP));
    K5_c = K2_c;

    dy(ind.Mp)  = K4_c * state(ind.AMp) + K1_c * SMC(flu.M) - (K2_c + K3_c) * state(ind.Mp);
    dy(ind.AMp) = K3_c * state(ind.Mp) + K6_c * state(ind.AM) - (K7_c + K5_c) * state(ind.AMp);
    dy(ind.AM)  = K5_c * state(ind.AMp) - (K7_c + K6_c) * state(ind.AM);

    %% ----- Radius dynamics -----
    F_r   = state(ind.AMp) + state(ind.AM);
    E_r   = Epas_r + F_r * (Eact_r - Epas_r);
    R0_r  = R0pas_r + F_r * R0pas_r * (0.6 - 1);
    dy(ind.R) = R0pas_r/nu_r * ( state(ind.R)*P_r/SMC(flu.h_r) - E_r * ((state(ind.R) - R0_r)/R0_r) );
end
