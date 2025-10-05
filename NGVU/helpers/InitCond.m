function STATES = InitCond()
%INITCOND Initial conditions for the merged NVC model.
% Returns a row vector STATES such that STATES(ind.<name>) is defined.
%
% Notes
% - The system is intended to be near steady state at t=0.
% - Vessel-dependent radius is selected from radii_list using global ii.

    % Globals / indices
    global ii radii_list;
    all_indices();   % defines structs 'ind' and 'flu'
    
    % ---- preallocate with last used index (35) ----
    NSTATES = 35;
    STATES  = zeros(1, NSTATES);

    % ---- Astrocyte/perivascular state (affects SMC via K_p / KIR) ----
    STATES(ind.R_k)      = 0.061e-6;     % m
    STATES(ind.N_Na_k)   = 0.99796e-3;   % uM*m
    STATES(ind.N_K_k)    = 5.52782e-3;   % uM*m
    STATES(ind.N_HCO3_k) = 0.58804e-3;   % uM*m
    STATES(ind.N_Cl_k)   = 0.32879e-3;   % uM*m
    STATES(ind.N_Na_s)   = 4.301041e-3;  % uM*m
    STATES(ind.N_K_s)    = 0.0807e-3;    % uM*m
    STATES(ind.N_HCO3_s) = 0.432552e-3;  % uM*m
    STATES(ind.K_p)      = 3e3;          % uM  (perivascular K+)
    STATES(ind.w_k)      = 0.1815e-3;    % [-] BK channel open prob.

    % ---- SMC ----
    STATES(ind.Ca_i)     = 0.1;          % uM
    STATES(ind.s_i)      = 0.1;          % uM
    STATES(ind.v_i)      = -60;          % mV
    STATES(ind.w_i)      = 0.1;          % [-]
    STATES(ind.I_i)      = 0.1;          % uM
    STATES(ind.K_i)      = 100e3;        % uM

    % ---- EC ----
    STATES(ind.Ca_j)     = 0.1;          % uM
    STATES(ind.s_j)      = 0.1;          % uM
    STATES(ind.v_j)      = -75;          % mV
    STATES(ind.I_j)      = 0.1;          % uM

    % ---- Crossbridge model ----
    STATES(ind.Mp)       = 0.25;
    STATES(ind.AMp)      = 0.25;
    STATES(ind.AM)       = 0.25;

    % ---- Vessel radius (per-vessel) ----
    if isempty(radii_list)
        error('InitCond: radii_list is empty. Set global radii_list in the driver before calling InitCond().');
    end
    if ii < 1 || ii > numel(radii_list)
        error('InitCond: ii=%d out of range for radii_list (size=%d).', ii, numel(radii_list));
    end
    STATES(ind.R) = radii_list(ii);

    % ---- Neuron/NO pathway initial conditions ----
    STATES(ind.Ca_n)        = 0.1;      % uM
    STATES(ind.nNOS_act_n)  = 0.29521;  % uM
    STATES(ind.NO_n)        = 0.0907;   % uM
    STATES(ind.NO_k)        = 0.0907;   % uM
    STATES(ind.O2_k)        = 200;      % uM
    STATES(ind.NO_i)        = 0.045;    % uM
    STATES(ind.NO_j)        = 0.045;    % uM
    STATES(ind.cGMP_i)      = 8;        % uM
    STATES(ind.eNOS_act_j)  = 0.7017;   % uM
    STATES(ind.E_b)         = 1;        % fraction
    STATES(ind.E_6c)        = 0;        % fraction
end
