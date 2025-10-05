% Quad-partite Synapse Model, Glutamate Neuron Dynamics
function [dydt] = GluN(~,y,param1,param2,param3,param4,param5, ...
                        param8,param9,param11,param14,param16, ...
                        Iapp,I_syn,Ggaba_syn,r_gabab,m_nmda,RRP, p_init)

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
gv = param5(2);
degG = param5(10);
v_glu = param8(1);
v_gaba = param8(2);
k_gaba = param8(3);
k_glu = param8(4);
tau_ip3 = param8(5);
np = param8(6);
tau_rec = param9(1);
tau_inact = param9(2);
F = param11(3);
g_nmda = param14(3);
V_nmda = param14(4);
Mg = param14(5);
g_gabab = param16(3);
V_GABAb = param16(4);

% Values set to equal input values
p = y(1);
m = y(2);
h = y(3);
n = y(4);
q = y(5);
mc= y(6);
v = y(7);
c = y(8);
cer=y(9);
R_syn=y(10);
E_syn=y(11);
G_syn=y(12);

% Activation/inactivation of ion-specific channels
an=0.01*((-v-60)/(exp((-v-60)/10)-1));
bn=0.125*exp((-v-70)/80);                                % K+
am=0.1*((-v-45)/(exp((-v-45)/10)-1));
bm=4*exp((-v-70)/18);                                    % Na+
ah=0.07*exp((-v-70)/20);
bh=1/(exp((-v-40)/10)+1);                                % Na+
aq=a2*d2*(p+d1)/(p+d3);
bq=a2*c;                                                 % IP3-R gating
mcinf=1/(1+exp((v_half-v)/kca));                         % VGCC channel

% Ion currents
ina=gna*(v-vna);
ik=gk*(v-vk);
il=gl*(v-vl);
ica=gc*(mc^2)*(v-vca);
I_pump=Ip*(c^2/(c^2+k_pump^2));
I_leak=v_leak*(v-vca);
I_gabab=g_gabab*r_gabab*(v-V_GABAb);
Mgv_glu =(1+exp((-0.062*v*Mg)/3.57))^(-1);      % Mg2+ ion binding dynamics
I_glu=g_nmda*m_nmda*(v-V_nmda)*Mgv_glu;         % Excitatory NMDA-R current (uA);

% IP3R Kinetics
minf=p/(p+d1);
ninf=c/(c+d5);                                     % IP3 & Ca2+ closing of IP3-R
jchan=c1*v1*(minf^3)*(ninf^3)*(q^3)*(c-cer);       % Ca++ flux through the IP3-R channel
jpump=v3*(c^2)/(k3^2+c^2);                         % SERCA pump flux
jleak=c1*v2*(c-cer);                               % Ca++ leak from the ER 
p_glu=v_glu*((G_syn^np)/(k_glu^np + G_syn^np));    % Extra-synaptic Glu-dependent IP3 production

% The Glutamate Neuron Dynamics
dydt = [(((p_init-p)*(tau_ip3))+p_glu+v_gaba*((Ggaba_syn^(np-.1)/(k_gaba^(np-.1) ...
        +Ggaba_syn^(np-0.1)))));                   % [IP3] in the nerve terminal
        (am*(1-m)-bm*m);                           % Na+ channel activation
        (ah*(1-h)-bh*h);                           % Na+ channel inactivation
        (an*(1-n)-bn*n);                           % K+ channel activation
        (aq*(1-q)-bq*q);                           % IP3-R gating variable
        (mcinf-mc)/taumc;                          % VGCC gating variable
    (Iapp-((m^3)*h*ina+(n^4)*ik+il) ...
        +I_gabab+I_glu);                          % Membrane potential (mV)
    (((-ica-I_pump-I_leak)*s_neur)/(2*F*v_neur) ...
        -jpump-jleak-jchan);                       % Intracellular [Ca++]i (nM)
    ((jchan+jleak+jpump)/c1);                      % [Ca++]er (nM)
    (((I_syn)/tau_rec)-(RRP*R_syn));             % Releasable Glu
    ((RRP*R_syn)-(E_syn/tau_inact));             % Effective Glu 
    (nv*gv*E_syn-degG*(G_syn))];                   % Glu in the cleft
end