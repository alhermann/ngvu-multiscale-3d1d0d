%% Merged Quadripartite + NVC driver
% Time is outermost → within each exchange window, loop vessels (1..50)
% For each vessel & window:
%   1) Advance quadripartite model up to the exchange time using RK4
%   2) Take ca, G_syn at the exchange instant
%   3) Run NVC ODE (ode15s) over that window with ca/G_syn held *constant*
%
% Assumes the following functions already exist on the MATLAB path:
% Markov, GluN, GabaN, Astrocyte, PostSynN, all_indices, all_constants,
% InitCond, quad_step_once, NVC_chunk, progress, all_fluxes, etc.

%% ----------------- Global configuration -----------------
% Global duration and exchange scheduling
T_TOTAL_SEC    = 60;             % total simulation time [s]
EXCH_POINTS    = 1200;           % number of exchange instants over 60 s (e.g., every 0.05 s)
% (You can change EXCH_POINTS to control exchange frequency.)

% Quadripartite integrator params (your original)
dt_ms          = 0.01;           % time-step for quadripartite RK4 [ms]

% NVC solver options (from your code)
RelTol         = 1e-3;
AbsTol         = 1e-3;
MaxStep        = 1;

global CASE ii radii_list
CASE           = 2;              % as in your NVC code

% ----------------- Vessel set -----------------
radii_list = [4.5e-06;3.5e-06;3.5e-06;3.5e-06;4.0e-06;4.0e-06;4.0e-06;3.5e-06; ...
              3.0e-06;3.0e-06;3.0e-06;3.0e-06;3.0e-06;2.5e-06;2.5e-06;2.5e-06; ...
              2.5e-06;2.5e-06;2.5e-06;2.5e-06;2.5e-06;2.5e-06;2.5e-06;2.5e-06; ...
              2.5e-06;2.5e-06;2.5e-06;2.5e-06;2.5e-06;2.5e-06;2.5e-06;2.0e-06; ...
              2.0e-06;2.0e-06;2.5e-06;2.5e-06;2.5e-06;2.5e-06;2.5e-06;2.5e-06; ...
              2.5e-06;3.5e-06;3.0e-06;2.5e-06;2.5e-06;2.5e-06;2.5e-06;2.5e-06; ...
              2.5e-06;2.5e-06];

% We keep independent state for quadripartite and NVC for each vessel.
VCOUNT = numel(radii_list);   % keep loops consistent with the radii you provide

% Output to microcirculation build dir (create if missing)
MIC_OUTDIR = fullfile('Microcirculation','Code_Upload','DUNE','dune-angiogenesis','build-cmake','src');
if ~isfolder(MIC_OUTDIR), mkdir(MIC_OUTDIR); end

% (Optional) legacy output folder not used for radius logging anymore
outdir = 'data';
if ~isfolder(outdir), mkdir(outdir); end

%% ----------------- Precompute exchange schedule -----------------
% Exchange times in seconds, strictly increasing 0 = t0 < t1 < ... <= tN
t_exch = linspace(0, T_TOTAL_SEC, EXCH_POINTS+1);  % length EXCH_POINTS+1

%% ----------------- Index and constants (shared) -----------------
% These define 'ind', 'flu', constants, etc. for both models
all_indices();
all_constants();

%% ----------------- Per-vessel initializations -----------------

% === Quadripartite state containers (per vessel) ===
quad = struct();
for ii = 1:VCOUNT
    rng(ii); % seed as in your original code
    
    % Geometry & electrical parameters (copied from your code)
    rad_neur=0.3141e-4; v_neur=(10^(-3))*(4/3)*pi*(rad_neur^3); s_neur=4*pi*(rad_neur^2);
    gna=120; gk=36; gl=0.3; p_ca=2.3e-9; rho_ca=(4*8.0659e+007); gc=p_ca*rho_ca;
    vl=-59.4; vna=45; vk=-82; vca=125; vr=-59.94;
    c1=0.185; v1=30e-3; v2=0.2374e-3; v3=90; k3=0.1e3; Ip=0.4; k_pump=100; v_leak=0.001022664392140;
    d1=0.13e3; d2=1.049e3; d3=0.9434e3; d5=0.08234e3; a2=0.2e-6;
    v_half=-17; kca=8.4; taumc=10;
    nv=2; gv=60; gabav=20; kon=0.3e-3; koff=3; kon_g=0.1e-2; koff_g=9; gamma=30; delta=8; degG=3; degGaba=1.5;
    sa1=50e3; sa2=5e3; sa3=0.85;
    sa_gaba1=50e3; sa_gaba2=5e3; sa_gaba3=0.77;
    v_glu=0.062; v_gaba=0.047; k_gaba=0.63e-3; k_glu=0.78e-3; tau_ip3=0.14e-3; np=0.7;
    tau_rec=800; tau_inact=3; tdr_glu=1000; tseq_glu=6.34;
    tau_rec_gaba=400; tau_inact_gaba=2; tdr_gaba=0; tseq_gaba=7.2;
    R=8.314; T=273.15+21; F=96487;
    R_in=0.7985e8; tau_mem=50; rad_ps=0.6e-4; vspine=(10^(-3))*(4/3)*pi*(rad_ps^3); sspine=4*pi*(rad_ps^2); v_post=-59.94; ks=100e-3;
    alpha_ampa=1.1; beta_ampa=0.67; g_ampa = 0.35e-6; V_ampa=0;
    alpha_nmda=0.72; beta_nmda=0.54; g_nmda=0.26e-6; V_nmda=0; Mg = 1;
    alpha_gaba=0.53; beta_gaba=0.18; g_gaba=0.745e-6; V_GABAa=-59.94;
    alpha_gabab=16; beta_gabab=6; g_gabab=0.445e-6; V_GABAb=-59.94;
    aNa=20; ac0=2e3; ac1=0.185; av1=6e-3; av2=0.11e-3; av3=0.9; ak3=0.1e3; ad1=0.13e3; ad2=1.049e3; ad3=0.9434e3; ad5=0.08234e3; aa2=0.2e-6;
    av_plcd=0.05; avgaba_plcd=0.37; ak_plcd=1.5e3; aK_plcd=0.1e3; ar5p=0.05e-3;
    av_plcb=0.5; avgaba_plcb=0.21; aK_R=1.3e-3; aK_R_gaba=0.95e-3; aK_P=10e-3; aK_P_gaba=10e-3;
    aK_pi=0.6e3; aK_pi_gaba=0.8e3; av_3K=2; aK_D=0.7e3; aK3=1e3;
    nva=12; nva_gaba=12; gva=20; gva_gaba=7; ak1=3.75e-6; ak2=2.5e-6; akk3=1.25e-5;
    ak_1=4e-4; ak_2=1e-3; ak_3=10e-3;
    adegG=10; atau_rec=800; atau_inact=3;
    adeg_gaba=8; atau_rec_gaba=1e3; atau_inact_gaba=2;

    % Bundle params in vectors used by your RHS calls
    param1 =[rad_neur,v_neur,s_neur];
    param2 = [gna,gk,gl,p_ca,rho_ca,gc,vl,vna,vk,vca,vr];
    param3 =[c1,v1,v2,v3,k3,Ip,k_pump,v_leak,d1,d2,d3,d5,a2];
    param4 = [v_half,kca,taumc];
    param5 =[nv,gv,gabav,kon,koff,kon_g,koff_g,gamma,delta,degG,degGaba];
    param6 =[sa1,sa2,sa3];
    param7 =[sa_gaba1,sa_gaba2,sa_gaba3];
    param8 =[v_glu,v_gaba,k_gaba,k_glu,tau_ip3,np];
    param9 =[tau_rec,tau_inact,tdr_glu,tseq_glu];
    param10 =[tau_rec_gaba,tau_inact_gaba,tdr_gaba,tseq_gaba];
    param11 =[R,T,F];
    param12 =[R_in,tau_mem,rad_ps,vspine,sspine,v_post,ks];
    param13 =[alpha_ampa,beta_ampa,g_ampa,V_ampa];
    param14 =[alpha_nmda,beta_nmda,g_nmda,V_nmda,Mg];
    param15 =[alpha_gaba,beta_gaba,g_gaba,V_GABAa];
    param16 =[alpha_gabab,beta_gabab,g_gabab,V_GABAb];
    param17 =[aNa,ac0,ac1,av1,av2,av3,ak3,ad1,ad2,ad3,ad5,aa2];
    param18 =[av_plcd,ak_plcd,aK_plcd,ar5p,av_plcb, ...
              avgaba_plcb,aK_R,aK_R_gaba,aK_P,aK_P_gaba,aK_pi, ...
              aK_pi_gaba,av_3K,aK_D,aK3];
    param19 =[nva,nva_gaba,gva,gva_gaba,ak1,ak2,akk3,ak_1,ak_2,ak_3];
    param20 =[adegG,atau_rec,atau_inact];
    param21 =[adeg_gaba,atau_rec_gaba,atau_inact_gaba];

    % Time bookkeeping for quadripartite model (internal ms clock)
    quad(ii).t_ms      = 0;        % current quadripartite 'time' in ms
    quad(ii).dt_ms     = dt_ms;
    quad(ii).params = struct('p1',param1,'p2',param2,'p3',param3,'p4',param4,'p5',param5, ...
                             'p6',param6,'p7',param7,'p8',param8,'p9',param9,'p10',param10, ...
                             'p11',param11,'p12',param12,'p13',param13,'p14',param14,'p15',param15, ...
                             'p16',param16,'p17',param17,'p18',param18,'p19',param19,'p20',param20,'p21',param21);

    % Initialize series (one-element, we will append per RK step)
    quad(ii).Iapp = 0; quad(ii).Iapp_gaba = 0;

    % Pre-synaptic (Glu)
    quad(ii).x1=Markov([1 0 0 0 0 0]); quad(ii).x2=Markov([1 0 0 0 0 0]);
    quad(ii).G_syn=1e-3; quad(ii).c=100; quad(ii).cer=1.30e6; quad(ii).v=-59.94;
    quad(ii).m=0.1; quad(ii).h=0.6; quad(ii).n=0.3; quad(ii).p=160; quad(ii).q=0.22; quad(ii).mc=0;
    quad(ii).R_syn=1; quad(ii).E_syn=0; quad(ii).I_syn=0; quad(ii).RRP=0; quad(ii).I_gabab=0;

    % Pre-synaptic (GABA)
    quad(ii).x1_gaba=Markov([1 0 0 0 0 0]); quad(ii).x2_gaba=Markov([1 0 0 0 0 0]);
    quad(ii).Ggaba_syn=1e-3; quad(ii).c_gaba=100; quad(ii).cer_gaba=1.30e6; quad(ii).V_gaba=-59.94;
    quad(ii).m_gaba=0.1; quad(ii).h_gaba=0.6; quad(ii).n_gaba=0.3; quad(ii).p_gaba=160; quad(ii).q_gaba=0.22; quad(ii).mc_gaba=0;
    quad(ii).Rgaba_syn=1; quad(ii).Egaba_syn=0; quad(ii).Igaba_syn=0; quad(ii).RRPgaba=0;

    % Astrocyte
    quad(ii).ca=100; quad(ii).ax=0.5; quad(ii).ap=160; quad(ii).aO1=0.48; quad(ii).aO2=0.2; quad(ii).aO3=5e-4;
    quad(ii).aE_syn=0; quad(ii).aI_syn=0; quad(ii).aR_syn=1; quad(ii).aG_syn=1e-3;
    quad(ii).aGaba_syn=1e-3; quad(ii).aEgaba_syn=0; quad(ii).aIgaba_syn=0; quad(ii).aRgaba_syn=1;

    % Post-synaptic
    quad(ii).V_post=-59.94; quad(ii).m_ampa=0; quad(ii).m_nmda=0; quad(ii).r_gaba=0; quad(ii).r_gabab=0;
    quad(ii).I_postampa=0; quad(ii).I_postnmda=0; quad(ii).I_postgaba=0; quad(ii).I_postsum=0;

    % carry refractory trackers
    quad(ii).tdr_glu = tdr_glu;
    quad(ii).tseq_glu = tseq_glu;
    quad(ii).tdr_gaba = tdr_gaba;
    quad(ii).tseq_gaba = tseq_gaba;

    quad(ii).p_init      = quad(ii).p;
    quad(ii).p_gaba_init = quad(ii).p_gaba;
end

%% ----------------- NVC per-vessel initialization -----------------
% Initialize NVC states and prepare per-vessel CSV files (cleared/created)
nvc = struct();
for ii = 1:VCOUNT
    state0 = InitCond();         % user-provided function
    nvc(ii).state  = state0(:);
    nvc(ii).t_last = 0;          % last absolute time [s] we integrated to

    % Prepare/clear per-vessel CSV in the requested folder
    vessel_csv = fullfile(MIC_OUTDIR, sprintf('NVU_Vessel_%d.csv', ii));
    if isfile(vessel_csv)
        delete(vessel_csv);      % clear previous runs
    else
        % touch the file so writematrix('append') always works
        writematrix([], vessel_csv);
    end
end

%% ----------------- Main loop: time outermost -----------------
fprintf('[MERGED] Running with %d exchange points over %g s (Δt ≈ %.5g s)\n', ...
        EXCH_POINTS, T_TOTAL_SEC, t_exch(2)-t_exch(1));

for k = 1:EXCH_POINTS
    t_a = t_exch(k);        % start time of window [s]
    t_b = t_exch(k+1);      % end time of window [s]

    % --- 1) Advance quadripartite for each vessel up to time t_b ---
    for ii = 1:VCOUNT
        t_target_ms = t_b * 1000;          % desired absolute quad time (ms)
        dt = quad(ii).dt_ms;
        steps = max(0, ceil((t_target_ms - quad(ii).t_ms) / dt));

        for s = 1:steps
            quad(ii) = quad_step_once(quad(ii));
        end

        % Record that we’ve advanced the internal ms clock
        quad(ii).t_ms = t_target_ms;
    end

    % --- 2) Run NVC over [t_a, t_b] with constant ca, G_syn per-vessel ---
    for ii = 1:VCOUNT
        % Take ca and G_syn at exchange instant (latest entries)
        ca_const   = quad(ii).ca;        % [nM]
        Gsyn_const = quad(ii).G_syn;     % [mM]

        % RHS for this window with constants (no time-indexing)
        ode_chunk = @(time,state) NVC_chunk(time, state, ca_const, Gsyn_const, T_TOTAL_SEC, t_a, t_b);

        options = odeset('OutputFcn',@progress,'Stats','on', ...
                         'RelTol',RelTol,'AbsTol',AbsTol,'MaxStep',MaxStep);

        % Integrate from current absolute time to next exchange time
        tspan = [nvc(ii).t_last, t_b];
        [t_chunk, state_chunk] = ode15s(ode_chunk, tspan, nvc(ii).state, options);

        % Append ONE row after each exchange:
        %   time [s], radius [micrometre]
        vessel_csv = fullfile(MIC_OUTDIR, sprintf('NVU_Vessel_%d.csv', ii));
        radius_um  = 1e6 * state_chunk(end, ind.R);
        writematrix([t_b, radius_um], vessel_csv, 'WriteMode','append');

        % Update last state/time
        nvc(ii).state  = state_chunk(end, :)';
        nvc(ii).t_last = t_b;
    end

    fprintf('[%4d/%4d] t = %.3f s → %.3f s  | completed all vessels\n', k, EXCH_POINTS, t_a, t_b);
end

fprintf('[DONE] Merged simulation finished for all %d vessels.\n', VCOUNT);
