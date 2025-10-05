% Quad-partite Synapse Model, GABA-Neuron Dynamics
function [dydt] = GabaN(~,y,param1,param2,param3,param4,param5, ...
                        param8,param10,param11,Iapp_gaba,Igaba_syn,RRPgaba,p_gaba_init)

% % Constants and parameters
v_neur=param1(2);
s_neur=param1(3);
gna = param2(1);
gk  = param2(2);
gl =  param2(3);
gc = param2(6);
vl = param2(7);
vna = param2(8);
vk = param2(9);
vca = param2(10);
c1 =param3(1);
v1 =param3(2);
v2 =param3(3);
v3 =param3(4);
k3 =param3(5);
Ip = param3(6);
k_pump = param3(7);
v_leak = param3(8);
d1 = param3(9);
d2 = param3(10);
d3 = param3(11);
d5 = param3(12);
a2 = param3(13);
v_half = param4(1);
kca = param4(2);
taumc = param4(3);
nv = param5(1);
gabav = param5(3);
degGaba = param5(11);
v_gaba = param8(2);
k_gaba = param8(3);
tau_ip3 = param8(5);
np = param8(6);
tau_rec_gaba = param10(1);
tau_inact_gaba = param10(2);
F = param11(3);

% Values set to equal input values
p_gaba = y(1);
m_gaba = y(2);
h_gaba = y(3);
n_gaba = y(4);
q_gaba = y(5);
mc_gaba= y(6);
V_gaba = y(7);
c_gaba = y(8);
cer_gaba=y(9);
Rgaba_syn=y(10);
Egaba_syn=y(11);
Ggaba_syn=y(12);

% Activation/inactivation variables
an_gaba=0.01*((-V_gaba-60)/(exp((-V_gaba-60)/10)-1));
bn_gaba=0.125*exp((-V_gaba-70)/80);
am_gaba=0.1*((-V_gaba-45)/(exp((-V_gaba-45)/10)-1));
bm_gaba=4*exp((-V_gaba-70)/18);
ah_gaba=0.07*exp((-V_gaba-70)/20);
bh_gaba=1/(exp((-V_gaba-40)/10)+1);
aq_gaba=a2*d2*(p_gaba+d1)/(p_gaba+d3);
bq_gaba=a2*c_gaba;                                          % IP3-R
mcinf_gaba=1/(1+exp((v_half-V_gaba)/kca));                  % Activity of the VGCC channel

% Ionic currents
ina_gaba=gna*(V_gaba-vna);
ik_gaba=gk*(V_gaba-vk);
il_gaba=gl*(V_gaba-vl);
ica_gaba=gc*(mc_gaba^2)*(V_gaba-vca);
I_pump_gaba=Ip*(c_gaba^2/(c_gaba^2+k_pump^2));
I_leakgaba=v_leak*(V_gaba-vca);

% IP3R Kinetics
minf_gaba=p_gaba/(p_gaba+d1);                               % Steady-state functions
ninf_gaba=c_gaba/(c_gaba+d5);                               % Steady-state functions
jchan_gaba=c1*v1*(minf_gaba^3)*(ninf_gaba^3)*....           % Ca++ flux from the ER to the cytosol
    (q_gaba^3)*(c_gaba-cer_gaba);
jpump_gaba=v3*(c_gaba^2)/(k3^2+c_gaba^2);                   % SERCA pump
jleak_gaba=c1*v2*(c_gaba-cer_gaba);                         % Ca++ leak from the ER into the cytosol

% The GABA-Neuron Dynamics
dydt = [(((p_gaba_init-p_gaba)*(tau_ip3))+ ...
    v_gaba*((Ggaba_syn^(np-.1)/(k_gaba^(np-.1) ...
    +Ggaba_syn^(np-0.1)))));                                % [IP3] produced
    (am_gaba*(1-m_gaba)-bm_gaba*m_gaba);                    % Na+ channel activation
    (ah_gaba*(1-h_gaba)-bh_gaba*h_gaba);                    % Na+ channel inactivation
    (an_gaba*(1-n_gaba)-bn_gaba*n_gaba);                    % K+ channel activation
    (aq_gaba*(1-q_gaba)-bq_gaba*q_gaba);                    % IP3-R gating variable
    (mcinf_gaba-mc_gaba)/taumc;                             % VGCC gating variable
    (Iapp_gaba-((m_gaba^3)*h_gaba*ina_gaba+ ...
    (n_gaba^4)*ik_gaba+il_gaba));                           % Membrane potential (mV)
    (((-ica_gaba-I_pump_gaba-I_leakgaba)*s_neur)/(2*F*v_neur) ...
    -jpump_gaba-jleak_gaba-jchan_gaba);                     % [Ca++]i (nM)
    ((jchan_gaba+jleak_gaba+jpump_gaba)/c1);                % [Ca++]er (nM)
    (((Igaba_syn)/tau_rec_gaba)-((RRPgaba)*Rgaba_syn));     % Fraction of releasable Glu
    (((RRPgaba)*Rgaba_syn)-(Egaba_syn/tau_inact_gaba));     % Fraction of effective Glu
    (nv*gabav*Egaba_syn-degGaba*(Ggaba_syn))];              % Glu in the cleft
end