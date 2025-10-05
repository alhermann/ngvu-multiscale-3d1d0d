function [NE, AC, SMC, EC] = all_fluxes(t, state, t_start, t_end, Glu_sc)
%ALL_FLUXES Algebraic relations/fluxes for NE, AC, SMC, EC.
% Uses: all_indices(), all_constants().
%
% Inputs:
%   t         : time [s]
%   state     : state vector (see ind.*)
%   t_start   : stimulus start [s]
%   t_end     : stimulus end [s]
%   Glu_sc    : glutamate in synaptic cleft [uM]
%
% Outputs:
%   NE, AC, SMC, EC : numeric arrays addressed with flu.* pointers

    all_indices();
    all_constants();

    % ========================= Neuron (NE) =========================
    NE = zeros(1, flu.tau_kn);  % preallocate up to last NE pointer

    % NMDA subunits open fractions (A/B)
    NE(flu.w_NR2_A) = Glu_sc / (K_m_A + Glu_sc);
    NE(flu.w_NR2_B) = Glu_sc / (K_m_B + Glu_sc);

    % Inward Ca2+ current per open NMDA receptor [fA]
    % (Relies on constants: v_n, G_M, P_Ca_P_M, Ca_ex, M, alpha_v, beta_v, R_gas, Temp, Farad)
    NE(flu.I_Ca)     = (4 * v_n * G_M * P_Ca_P_M * Ca_ex / M) / (1 + exp(alpha_v * (v_n + beta_v))) ...
                     * (exp(2 * v_n * Farad / (R_gas * Temp))) / (1 - exp(2 * v_n * Farad / (R_gas * Temp)));

    % Total Ca2+ current across all open NMDA receptors
    NE(flu.I_Ca_tot) = (n_NR2_A * NE(flu.w_NR2_A) + n_NR2_B * NE(flu.w_NR2_B)) * NE(flu.I_Ca);

    % Calmodulin binding polynomials
    NE(flu.phi_mc)   = 1 + Q_1 * state(ind.Ca_n) + Q_1*Q_2 * state(ind.Ca_n)^2 ...
                        + Q_1*Q_2*Q_3 * state(ind.Ca_n)^3 + Q_1*Q_2*Q_3*Q_4 * state(ind.Ca_n)^4;

    % Expected Ca bound per CaM
    NE(flu.m_c)      = (state(ind.Ca_n)) / NE(flu.phi_mc) ...
                      * (Q_1 + 2*Q_1*Q_2*state(ind.Ca_n) + 3*Q_1*Q_2*Q_3*state(ind.Ca_n)^2 ...
                      + 4*Q_1*Q_2*Q_3*Q_4*state(ind.Ca_n)^3);

    % NO kinetics (neuron)
    NE(flu.p_NO_n)   = V_max_NO_n * state(ind.nNOS_act_n) * O2_n / (K_m_O2_n + O2_n) * LArg_n / (K_m_LArg_n + LArg_n);
    NE(flu.c_NO_n)   = k_O2_n * state(ind.NO_n)^2 * O2_n;

    % Diffusive time and flux (neuron ↔ astrocyte)
    NE(flu.tau_kn)   = x_nk^2 / (2 * D_c_NO);
    NE(flu.d_NO_n)   = (state(ind.NO_k) - state(ind.NO_n)) / NE(flu.tau_kn);

    % ======================= Astrocyte (AC) ========================
    AC = zeros(1, flu.tau_ki);  % preallocate up to last AC pointer

    % Geometry
    AC(flu.R_s)    = R_tot - state(ind.R_k);               % [m]
    AC(flu.N_Cl_s) = state(ind.N_Na_s) + state(ind.N_K_s) - state(ind.N_HCO3_s);

    % Convert extensive amounts (uM*m) to concentrations (uM) with compartment widths
    AC(flu.Na_k)   = negCheck(state(ind.N_Na_k),   state(ind.R_k));
    AC(flu.K_k)    = negCheck(state(ind.N_K_k),    state(ind.R_k));
    AC(flu.HCO3_k) = negCheck(state(ind.N_HCO3_k), state(ind.R_k));
    AC(flu.Cl_k)   = negCheck(state(ind.N_Cl_k),   state(ind.R_k));
    AC(flu.Na_s)   = negCheck(state(ind.N_Na_s),   AC(flu.R_s));
    AC(flu.K_s)    = negCheck(state(ind.N_K_s),    AC(flu.R_s));
    AC(flu.HCO3_s) = negCheck(state(ind.N_HCO3_s), AC(flu.R_s));
    AC(flu.Cl_s)   = negCheck(AC(flu.N_Cl_s),      AC(flu.R_s));

    % Nernst potentials (V)
    AC(flu.E_Na_k)  = (R_gas * Temp) / (z_Na * Farad) * log(AC(flu.Na_s)/AC(flu.Na_k));
    AC(flu.E_K_k)   = (R_gas * Temp) / (z_K  * Farad) * log(AC(flu.K_s) /AC(flu.K_k ));
    AC(flu.E_Cl_k)  = (R_gas * Temp) / (z_Cl * Farad) * log(AC(flu.Cl_s)/AC(flu.Cl_k));
    AC(flu.E_NBC_k) = (R_gas * Temp) / (z_NBC* Farad) ...
                    * log((AC(flu.Na_s)*AC(flu.HCO3_s)^2)/(AC(flu.Na_k)*AC(flu.HCO3_k)^2));
    AC(flu.E_BK_k)  = (R_gas * Temp) / (z_K  * Farad) * log(state(ind.K_p)/AC(flu.K_k));

    % Na/K pump
    AC(flu.J_NaK_k) = J_NaK_max * Hill(AC(flu.Na_k), K_Na_k, 1.5) * Hill(AC(flu.K_s), K_K_s, 1);

    % Membrane potential (algebraic balance)
    AC(flu.v_k) = ( g_Na_k  * AC(flu.E_Na_k) ...
                  + g_K_k   * AC(flu.E_K_k ) ...
                  + g_Cl_k  * AC(flu.E_Cl_k) ...
                  + g_NBC_k * AC(flu.E_NBC_k) ...
                  - AC(flu.J_NaK_k)*Farad/unitcon ...
                  + g_BK_k  * state(ind.w_k) * AC(flu.E_BK_k) ) ...
                  /(g_Na_k + g_K_k + g_Cl_k + g_NBC_k + g_BK_k*state(ind.w_k));

    % Transporters/channels (uM*m s^-1)
    AC(flu.J_KCC1_k ) = 0.0;

    AC(flu.J_NBC_k  ) = g_NBC_k  / Farad * (AC(flu.v_k) - AC(flu.E_NBC_k)) * unitcon;
    AC(flu.J_NKCC1_k) = 0.0;
    AC(flu.J_Na_k   ) = g_Na_k   / Farad * (AC(flu.v_k) - AC(flu.E_Na_k)) * unitcon;
    AC(flu.J_K_k    ) = g_K_k    / Farad * (AC(flu.v_k) - AC(flu.E_K_k )) * unitcon;
    AC(flu.J_BK_k   ) = g_BK_k   / Farad * state(ind.w_k) * (AC(flu.v_k) - AC(flu.E_BK_k)) * unitcon;

    % BK gating
    AC(flu.w_inf) = 0.5 * (1 + tanh((AC(flu.v_k) + v_6) / (v_4)));
    AC(flu.phi_w) = psi_w * cosh((AC(flu.v_k) + v_6) / (2*v_4));

    % NO in astrocyte
    AC(flu.p_NO_k) = 0;
    AC(flu.c_NO_k) = k_O2_k * state(ind.NO_k)^2 * state(ind.O2_k);

    % Diffusive coupling (NE↔AC and AC↔SMC)
    AC(flu.tau_ki) = x_ki^2 / (2 * D_c_NO);
    AC(flu.d_NO_k) = (state(ind.NO_n) - state(ind.NO_k)) / NE(flu.tau_kn) ...
                   + (state(ind.NO_i) - state(ind.NO_k)) / AC(flu.tau_ki);

    % ========================== SMC ==========================
    SMC = zeros(1, flu.tau_ij);

    SMC(flu.M)         = 1 - state(ind.Mp) - state(ind.AM) - state(ind.AMp);
    SMC(flu.E_K_i)     = (R_gas * Temp) / (z_K * Farad) * unitcon * log(state(ind.K_p)/state(ind.K_i));
    SMC(flu.h_r)       = -state(ind.R) + sqrt(state(ind.R)^2 + 2*rb_r*h0_r + h0_r^2);

    % Electrical/chemical coupling with EC
    SMC(flu.v_coup_i)  = - g_hat   * (state(ind.v_i)   - state(ind.v_j));
    SMC(flu.Ca_coup_i) = - p_hat   * (state(ind.Ca_i)  - state(ind.Ca_j));
    SMC(flu.IP3_coup_i)= - p_hatIP3* (state(ind.I_i)   - state(ind.I_j));

    % Ca handling / currents
    SMC(flu.rho_i)         = 1;
    SMC(flu.J_IP3_i)       = Fmax_i   * (state(ind.I_i)^2) / (Kr_i^2 + state(ind.I_i)^2);
    SMC(flu.J_SRuptake_i)  = B_i      * (state(ind.Ca_i)^2) / (state(ind.Ca_i)^2 + cb_i^2);
    SMC(flu.J_CICR_i)      = C_i      * (state(ind.s_i)^2) / (sc_i^2 + state(ind.s_i)^2) ...
                                       * (state(ind.Ca_i)^4) / (cc_i^4 + state(ind.Ca_i)^4);
    SMC(flu.J_extrusion_i) = D_i      * state(ind.Ca_i) * ( 1 + (state(ind.v_i) - vd_i) / Rd_i );
    SMC(flu.J_leak_i)      = L_i      * state(ind.s_i);
    SMC(flu.J_VOCC_i)      = G_Ca     * (state(ind.v_i) - v_Ca1) / (1 + exp(-(state(ind.v_i) - v_Ca2) / R_Ca));
    SMC(flu.J_NaCa_i)      = G_NaCa   * state(ind.Ca_i) * (state(ind.v_i) - v_NaCa) / (state(ind.Ca_i) + c_NaCa);
    SMC(flu.J_NaK_i)       = F_NaK;
    SMC(flu.J_Cl_i)        = G_Cl     * (state(ind.v_i) - v_Cl);
    SMC(flu.J_K_i)         = G_K      * state(ind.w_i) * (state(ind.v_i) - vK_i);
    SMC(flu.J_degrad_i)    = k_i      * state(ind.I_i);
    SMC(flu.J_stretch_i)   = G_stretch/(1+exp(-alpha1*(P_str*state(ind.R)/SMC(flu.h_r) - sig0))) * (state(ind.v_i) - Esac);

    % KIR (depends on K_p!)
    SMC(flu.v_KIR_i)  = z_1 * state(ind.K_p)/unitcon + z_2; % mV
    SMC(flu.G_KIR_i)  = exp(z_5 * state(ind.v_i) + z_3 * state(ind.K_p)/unitcon + z_4);
    SMC(flu.J_KIR_i)  = F_il/gam * SMC(flu.G_KIR_i) * (state(ind.v_i) - SMC(flu.v_KIR_i));

    % NO / cGMP in SMC
    SMC(flu.p_NO_i)    = 0;
    SMC(flu.c_NO_i)    = k_dno * state(ind.NO_i);
    SMC(flu.tau_ij)    = (x_ij)^2 / (2 * D_c_NO);
    SMC(flu.d_NO_i)    = (state(ind.NO_k) - state(ind.NO_i)) / AC(flu.tau_ki) + (state(ind.NO_j) - state(ind.NO_i)) / SMC(flu.tau_ij);
    SMC(flu.k_4)       = C_4 * (state(ind.cGMP_i))^(m_4);
    SMC(flu.E_5c)      = 1 - state(ind.E_b) - state(ind.E_6c);
    SMC(flu.R_cGMP)    = (state(ind.cGMP_i)^2) / ((K_m_mlpc)^2 + (state(ind.cGMP_i)^2));
    SMC(flu.K_5c)      = SMC(flu.K_2c);  % alias used below
    SMC(flu.K_2c)      = delta_i * (k_mlpc_b + k_mlpc_c * SMC(flu.R_cGMP));
    SMC(flu.c_w_i)     = 1 / (epsilon_i + alpha_i * exp(gam_i * state(ind.cGMP_i)));
    SMC(flu.K_act_i)   = (state(ind.Ca_i) + SMC(flu.c_w_i)) / (state(ind.Ca_i) + (SMC(flu.c_w_i))^2 + beta_i * exp(v_Ca3_i - state(ind.v_i) / R_K_i));
    SMC(flu.V_max_pde) = k_pde * state(ind.cGMP_i);

    % Activation used also in w_i kinetic
    SMC(flu.Kactivation_i) = (state(ind.Ca_i) + SMC(flu.c_w_i))^2 / ( (state(ind.Ca_i) + SMC(flu.c_w_i))^2 + bet * exp(-(state(ind.v_i) - v_Ca3)/R_K) );

    % =========================== EC ===========================
    EC = zeros(1, flu.d_NO_j);

    % Coupling
    EC(flu.v_coup_j)   = - g_hat   * (state(ind.v_j)  - state(ind.v_i));
    EC(flu.Ca_coup_j)  = - p_hat   * (state(ind.Ca_j) - state(ind.Ca_i));
    EC(flu.IP3_coup_j) = - p_hatIP3* (state(ind.I_j)  - state(ind.I_i));

    % Ca handling / currents
    EC(flu.rho_j)        = 1;
    EC(flu.J_0_j)        = J0_j;
    EC(flu.J_IP3_j)      = Fmax_j * (state(ind.I_j)^2) / (Kr_j^2 + state(ind.I_j)^2);
    EC(flu.J_ERuptake_j) = B_j    * (state(ind.Ca_j)^2) / (state(ind.Ca_j)^2 + cb_j^2);
    EC(flu.J_CICR_j)     = C_j    * (state(ind.s_j)^2) / (sc_j^2 + state(ind.s_j)^2) ...
                                     * (state(ind.Ca_j)^4) / (cc_j^4 + state(ind.Ca_j)^4);
    EC(flu.J_extrusion_j)= D_j    * state(ind.Ca_j);
    EC(flu.J_leak_j)     = L_j    * state(ind.s_j);

    EC(flu.J_cation_j)   = G_cat * (E_Ca - state(ind.v_j)) * 0.5 * (1 + tanh((log10(state(ind.Ca_j)) - m3cat)/m4cat));
    EC(flu.J_BKCa_j)     = 0.4/2 * (1 + tanh( ((log10(state(ind.Ca_j)) - c)*(state(ind.v_j) - b) - a1) ...
                            / ( m3b*(state(ind.v_j) + a2*(log10(state(ind.Ca_j)) - c) - b)^2 + m4b )));
    EC(flu.J_SKCa_j)     = 0.6/2 * (1 + tanh( (log10(state(ind.Ca_j)) - m3s) / m4s ));
    EC(flu.J_K_j)        = G_tot * (state(ind.v_j) - vK_j) * (EC(flu.J_BKCa_j) + EC(flu.J_SKCa_j));
    EC(flu.J_R_j)        = G_R   * (state(ind.v_j) - v_rest);
    EC(flu.J_degrad_j)   = k_j   * state(ind.I_j);
    EC(flu.J_stretch_j)  = G_stretch/(1+exp(-alpha1*(P_str*state(ind.R)/SMC(flu.h_r) - sig0))) * (state(ind.v_j) - Esac);

    % Wall shear stress / NO in EC
    EC(flu.tau_wss) = state(ind.R) / 2 * DeltaPL;
    EC(flu.W_wss)   = W_0 * (EC(flu.tau_wss) + sqrt(16*(delta_wss)^2 + EC(flu.tau_wss)) - 4*delta_wss)^2 ...
                    / (EC(flu.tau_wss) + sqrt(16*(delta_wss)^2 - (EC(flu.tau_wss)^2)));
    EC(flu.F_wss)   = 1/(1 + alpha_wss*exp(-EC(flu.W_wss))) - 1/(1 + alpha_wss);
    EC(flu.p_NO_j)  = (V_max_NO_j * state(ind.eNOS_act_j) * O2_j) / (K_m_O2_j + O2_j) * LArg_j / (K_m_LArg_j + LArg_j);
    EC(flu.c_NO_j)  = k_O2_j * (state(ind.NO_j))^2 * O2_j;
    EC(flu.d_NO_j)  = (state(ind.NO_i) - state(ind.NO_j)) / SMC(flu.tau_ij) - (4 * D_c_NO * state(ind.NO_j)) / (state(ind.R))^2;
end
