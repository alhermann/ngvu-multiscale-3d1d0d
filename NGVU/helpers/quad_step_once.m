%% ===================== Helper: one RK4 step of quadripartite =====================
function Q = quad_step_once(Q)
% Executes one dt_ms step of your quadripartite loop for ONE vessel.
% This is a faithful transcription of your per-step math (Glu, GABA, Astro, Post).

% Unpack parameters for readability
p1=Q.params.p1; p2=Q.params.p2; p3=Q.params.p3; p4=Q.params.p4; p5=Q.params.p5;
p6=Q.params.p6; p7=Q.params.p7; p8=Q.params.p8; p9=Q.params.p9; p10=Q.params.p10; p11=Q.params.p11;
p12=Q.params.p12; p13=Q.params.p13; p14=Q.params.p14; p15=Q.params.p15; p16=Q.params.p16;
p17=Q.params.p17; p18=Q.params.p18; p19=Q.params.p19; p20=Q.params.p20; p21=Q.params.p21;

dt = Q.dt_ms;                  % ms
t  = Q.t_ms;                   % current absolute time (ms)

%% ===== Glutamate neuron =====
% Stimulus: 5 Hz, 4 ms width every 10,000 ms
if mod(t,10000)>=0 && mod(t,10000)<=4
    Iapp = 10;
else
    Iapp = 0;
end

% Spontaneous release rate (sa parameters)
lambda = p6(3)/(1+exp((p6(1)-Q.c)/p6(2)));
poisson_rand_lambda = poissrnd(lambda);

% Transition probabilities for Ca-sensor
kon  = p5(4);  koff  = p5(5);  gamma = p5(8);  delta = p5(9);
a21=5*kon*Q.c*dt;   a12=koff*dt;
a32=4*kon*Q.c*dt;   a23=2*koff*dt;
a43=3*kon*Q.c*dt;   a34=3*koff*dt;
a54=2*kon*Q.c*dt;   a45=4*koff*dt;
a65=kon*Q.c*dt;     a56=5*koff*dt;
a76=gamma*dt;       a67=delta*dt;
D0=1-a21; D1=1-a12-a32; D2=1-a23-a43; D3=1-a34-a54; D4=1-a45-a65; D5=1-a56-a76; D6=1-a67;
P=[D0,a12,0,0,0,0,0;
   a21,D1,a23,0,0,0,0;
   0,a32,D2,a34,0,0,0;
   0,0,a43,D3,a45,0,0;
   0,0,0,a54,D4,a56,0;
   0,0,0,0,a65,D5,a67;
   0,0,0,0,0,a76,D6];

x1_next = Markov(P(:,Q.x1));
x2_next = Markov(P(:,Q.x2));

% Refractory and RRP logic
vr = p2(11);
if poisson_rand_lambda==2 && abs(Q.tdr_glu - t)>Q.tseq_glu && Q.v<=vr
    RRP = 1;
elseif poisson_rand_lambda==1 && abs(Q.tdr_glu - t)>Q.tseq_glu && Q.v<=vr
    RRP = 0.5; Q.tdr_glu=t;
elseif (Q.x1==7 && Q.x2~=7 && abs(Q.tdr_glu - t)>Q.tseq_glu) && Q.v>vr
    RRP = 0.5; Q.tdr_glu=t;
elseif (Q.x2==7 && Q.x1~=7 && abs(Q.tdr_glu - t)>Q.tseq_glu) && Q.v>vr
    RRP = 0.5; Q.tdr_glu=t;
elseif (Q.x1==7 && Q.x2==7 && abs(Q.tdr_glu - t)>Q.tseq_glu) && Q.v>vr
    RRP = 1; Q.tdr_glu=t;
else
    RRP = 0;
end

% RK4 for Glu neuron
K1=dt*GluN(0,[Q.p;Q.m;Q.h;Q.n;Q.q;Q.mc;Q.v;Q.c;Q.cer;Q.R_syn;Q.E_syn;Q.G_syn], ...
            p1,p2,p3,p4,p5,p8,p9,p11,p14,p16,Iapp,Q.I_syn,Q.Ggaba_syn,Q.r_gabab,Q.m_nmda,RRP,Q.p_init);
K2=dt*GluN(0.5*dt,[Q.p+0.5*K1(1);Q.m+0.5*K1(2);Q.h+0.5*K1(3);Q.n+0.5*K1(4); ...
                   Q.q+0.5*K1(5);Q.mc+0.5*K1(6);Q.v+0.5*K1(7);Q.c+0.5*K1(8);Q.cer+0.5*K1(9); ...
                   Q.R_syn+0.5*K1(10);Q.E_syn+0.5*K1(11);Q.G_syn+0.5*K1(12)], ...
            p1,p2,p3,p4,p5,p8,p9,p11,p14,p16,Iapp,Q.I_syn,Q.Ggaba_syn,Q.r_gabab,Q.m_nmda,RRP,Q.p_init);
K3=dt*GluN(0.5*dt,[Q.p+0.5*K2(1);Q.m+0.5*K2(2);Q.h+0.5*K2(3);Q.n+0.5*K2(4); ...
                   Q.q+0.5*K2(5);Q.mc+0.5*K2(6);Q.v+0.5*K2(7);Q.c+0.5*K2(8);Q.cer+0.5*K2(9); ...
                   Q.R_syn+0.5*K2(10);Q.E_syn+0.5*K2(11);Q.G_syn+0.5*K2(12)], ...
            p1,p2,p3,p4,p5,p8,p9,p11,p14,p16,Iapp,Q.I_syn,Q.Ggaba_syn,Q.r_gabab,Q.m_nmda,RRP,Q.p_init);
K4=dt*GluN(dt,[Q.p+K3(1);Q.m+K3(2);Q.h+K3(3);Q.n+K3(4); ...
               Q.q+K3(5);Q.mc+K3(6);Q.v+K3(7);Q.c+K3(8);Q.cer+K3(9); ...
               Q.R_syn+K3(10);Q.E_syn+K3(11);Q.G_syn+K3(12)], ...
            p1,p2,p3,p4,p5,p8,p9,p11,p14,p16,Iapp,Q.I_syn,Q.Ggaba_syn,Q.r_gabab,Q.m_nmda,RRP,Q.p_init);

Q.p    = Q.p    + (K1(1)+2*K2(1)+2*K3(1)+K4(1))/6;
Q.m    = Q.m    + (K1(2)+2*K2(2)+2*K3(2)+K4(2))/6;
Q.h    = Q.h    + (K1(3)+2*K2(3)+2*K3(3)+K4(3))/6;
Q.n    = Q.n    + (K1(4)+2*K2(4)+2*K3(4)+K4(4))/6;
Q.q    = Q.q    + (K1(5)+2*K2(5)+2*K3(5)+K4(5))/6;
Q.mc   = Q.mc   + (K1(6)+2*K2(6)+2*K3(6)+K4(6))/6;
Q.v    = Q.v    + (K1(7)+2*K2(7)+2*K3(7)+K4(7))/6;
Q.c    = Q.c    + (K1(8)+2*K2(8)+2*K3(8)+K4(8))/6;
Q.cer  = Q.cer  + (K1(9)+2*K2(9)+2*K3(9)+K4(9))/6;
Q.R_syn= Q.R_syn+ (K1(10)+2*K2(10)+2*K3(10)+K4(10))/6;
Q.E_syn= Q.E_syn+ (K1(11)+2*K2(11)+2*K3(11)+K4(11))/6;
Q.I_syn= 1 - Q.R_syn - Q.E_syn;
Q.G_syn= Q.G_syn+ (K1(12)+2*K2(12)+2*K3(12)+K4(12))/6;

Q.Iapp = Iapp;  % store last applied current
Q.x1   = x1_next;
Q.x2   = x2_next;
Q.RRP  = RRP;

%% ===== GABA neuron =====
if mod(t,8000)>=0 && mod(t,8000)<=4
    Iapp_gaba = 10;
else
    Iapp_gaba = 0;
end

lambda_gaba = p7(3)/(1+exp((p7(1)-Q.c_gaba)/p7(2)));
poisson_rand_lambda = poissrnd(lambda_gaba);

kon_g = p5(6); koff_g = p5(7);
a21=5*kon_g*Q.c_gaba*dt; a12=koff_g*dt;
a32=4*kon_g*Q.c_gaba*dt; a23=2*koff_g*dt;
a43=3*kon_g*Q.c_gaba*dt; a34=3*koff_g*dt;
a54=2*kon_g*Q.c_gaba*dt; a45=4*koff_g*dt;
a65=kon_g*Q.c_gaba*dt;   a56=5*koff_g*dt;
a76=p5(8)*dt;            a67=p5(9)*dt;
D0=1-a21; D1=1-a12-a32; D2=1-a23-a43; D3=1-a34-a54; D4=1-a45-a65; D5=1-a56-a76; D6=1-a67;
P=[D0,a12,0,0,0,0,0;
   a21,D1,a23,0,0,0,0;
   0,a32,D2,a34,0,0,0;
   0,0,a43,D3,a45,0,0;
   0,0,0,a54,D4,a56,0;
   0,0,0,0,a65,D5,a67;
   0,0,0,0,0,a76,D6];

x1g_next = Markov(P(:,Q.x1_gaba));
x2g_next = Markov(P(:,Q.x2_gaba));

if poisson_rand_lambda==2 && abs(Q.tdr_gaba - t)>Q.tseq_gaba && Q.V_gaba<=vr
    RRPg = 1;
elseif poisson_rand_lambda==1 && abs(Q.tdr_gaba - t)>Q.tseq_gaba && Q.V_gaba<=vr
    RRPg = 0.5; Q.tdr_gaba=t;
elseif (Q.x1_gaba==7 && Q.x2_gaba~=7 && abs(Q.tdr_gaba - t)>Q.tseq_gaba) && Q.V_gaba>vr
    RRPg = 0.5; Q.tdr_gaba=t;
elseif (Q.x2_gaba==7 && Q.x1_gaba~=7 && abs(Q.tdr_gaba - t)>Q.tseq_gaba) && Q.V_gaba>vr
    RRPg = 0.5; Q.tdr_gaba=t;
elseif (Q.x1_gaba==7 && Q.x2_gaba==7 && abs(Q.tdr_gaba - t)>Q.tseq_gaba) && Q.V_gaba>vr
    RRPg = 1; Q.tdr_gaba=t;
else
    RRPg = 0;
end

K1=dt*GabaN(0,[Q.p_gaba;Q.m_gaba;Q.h_gaba;Q.n_gaba;Q.q_gaba; ...
               Q.mc_gaba;Q.V_gaba;Q.c_gaba;Q.cer_gaba;Q.Rgaba_syn;Q.Egaba_syn;Q.Ggaba_syn], ...
            p1,p2,p3,p4,p5,p8,p10,p11,Iapp_gaba,Q.Igaba_syn,RRPg,Q.p_gaba_init);
K2=dt*GabaN(0.5*dt,[Q.p_gaba+0.5*K1(1);Q.m_gaba+0.5*K1(2);Q.h_gaba+0.5*K1(3);Q.n_gaba+0.5*K1(4); ...
                    Q.q_gaba+0.5*K1(5);Q.mc_gaba+0.5*K1(6);Q.V_gaba+0.5*K1(7);Q.c_gaba+0.5*K1(8); ...
                    Q.cer_gaba+0.5*K1(9);Q.Rgaba_syn+0.5*K1(10);Q.Egaba_syn+0.5*K1(11);Q.Ggaba_syn+0.5*K1(12)], ...
            p1,p2,p3,p4,p5,p8,p10,p11,Iapp_gaba,Q.Igaba_syn,RRPg,Q.p_gaba_init);
K3=dt*GabaN(0.5*dt,[Q.p_gaba+0.5*K2(1);Q.m_gaba+0.5*K2(2);Q.h_gaba+0.5*K2(3);Q.n_gaba+0.5*K2(4); ...
                    Q.q_gaba+0.5*K2(5);Q.mc_gaba+0.5*K2(6);Q.V_gaba+0.5*K2(7);Q.c_gaba+0.5*K2(8); ...
                    Q.cer_gaba+0.5*K2(9);Q.Rgaba_syn+0.5*K2(10);Q.Egaba_syn+0.5*K2(11);Q.Ggaba_syn+0.5*K2(12)], ...
            p1,p2,p3,p4,p5,p8,p10,p11,Iapp_gaba,Q.Igaba_syn,RRPg,Q.p_gaba_init);
K4=dt*GabaN(dt,[Q.p_gaba+K3(1);Q.m_gaba+K3(2);Q.h_gaba+K3(3);Q.n_gaba+K3(4); ...
                Q.q_gaba+K3(5);Q.mc_gaba+K3(6);Q.V_gaba+K3(7);Q.c_gaba+K3(8); ...
                Q.cer_gaba+K3(9);Q.Rgaba_syn+K3(10);Q.Egaba_syn+K3(11);Q.Ggaba_syn+K3(12)], ...
            p1,p2,p3,p4,p5,p8,p10,p11,Iapp_gaba,Q.Igaba_syn,RRPg,Q.p_gaba_init);

Q.p_gaba   = Q.p_gaba   + (K1(1)+2*K2(1)+2*K3(1)+K4(1))/6;
Q.m_gaba   = Q.m_gaba   + (K1(2)+2*K2(2)+2*K3(2)+K4(2))/6;
Q.h_gaba   = Q.h_gaba   + (K1(3)+2*K2(3)+2*K3(3)+K4(3))/6;
Q.n_gaba   = Q.n_gaba   + (K1(4)+2*K2(4)+2*K3(4)+K4(4))/6;
Q.q_gaba   = Q.q_gaba   + (K1(5)+2*K2(5)+2*K3(5)+K4(5))/6;
Q.mc_gaba  = Q.mc_gaba  + (K1(6)+2*K2(6)+2*K3(6)+K4(6))/6;
Q.V_gaba   = Q.V_gaba   + (K1(7)+2*K2(7)+2*K3(7)+K4(7))/6;
Q.c_gaba   = Q.c_gaba   + (K1(8)+2*K2(8)+2*K3(8)+K4(8))/6;
Q.cer_gaba = Q.cer_gaba + (K1(9)+2*K2(9)+2*K3(9)+K4(9))/6;
Q.Rgaba_syn= Q.Rgaba_syn+ (K1(10)+2*K2(10)+2*K3(10)+K4(10))/6;
Q.Egaba_syn= Q.Egaba_syn+ (K1(11)+2*K2(11)+2*K3(11)+K4(11))/6;
Q.Igaba_syn= 1 - Q.Rgaba_syn - Q.Egaba_syn;
Q.Ggaba_syn= Q.Ggaba_syn+ (K1(12)+2*K2(12)+2*K3(12)+K4(12))/6;

Q.Iapp_gaba = Iapp_gaba;
Q.x1_gaba   = x1g_next;
Q.x2_gaba   = x2g_next;
Q.RRPgaba   = RRPg;

%% ===== Astrocyte =====
aa2 = p17(12); ad2=p17(10); ad1=p17(8); ad3=p17(11);
aaq=aa2*ad2*(Q.ap+ad1)/(Q.ap+ad3);
abq=aa2*Q.ca;
au1=rand; au2=rand;
aa=((aaq*(1-Q.ax))-abq*Q.ax)/p17(1);
dW=sqrt(abs(-(2*dt*(aa)*log(au1))))*cos(2*pi*au2);

if Q.ax>=0 && Q.ax<=1
    ax_next = dt*(aaq*(1-Q.ax)-abq*Q.ax)+Q.ax+dW;
else
    ax_next = dt*(aaq*(1-Q.ax)-abq*Q.ax)+Q.ax;
end

K1=dt*Astrocyte(0,[Q.ca;Q.ap;Q.aO1;Q.aO2;Q.aO3;Q.aR_syn;Q.aE_syn;Q.aRgaba_syn; ...
                   Q.aEgaba_syn;Q.aG_syn;Q.aGaba_syn],p8,p17,p18,p19,p20,p21,Q.aI_syn,Q.aIgaba_syn,Q.ax,Q.G_syn,Q.Ggaba_syn);
K2=dt*Astrocyte(0.5*dt,[Q.ca+0.5*K1(1);Q.ap+0.5*K1(2);Q.aO1+0.5*K1(3);Q.aO2+0.5*K1(4);Q.aO3+0.5*K1(5); ...
                        Q.aR_syn+0.5*K1(6);Q.aE_syn+0.5*K1(7);Q.aRgaba_syn+0.5*K1(8);Q.aEgaba_syn+0.5*K1(9); ...
                        Q.aG_syn+0.5*K1(10);Q.aGaba_syn+0.5*K1(11)],p8,p17,p18,p19,p20,p21,Q.aI_syn,Q.aIgaba_syn,Q.ax,Q.G_syn,Q.Ggaba_syn);
K3=dt*Astrocyte(0.5*dt,[Q.ca+0.5*K2(1);Q.ap+0.5*K2(2);Q.aO1+0.5*K2(3);Q.aO2+0.5*K2(4);Q.aO3+0.5*K2(5); ...
                        Q.aR_syn+0.5*K2(6);Q.aE_syn+0.5*K2(7);Q.aRgaba_syn+0.5*K2(8);Q.aEgaba_syn+0.5*K2(9); ...
                        Q.aG_syn+0.5*K2(10);Q.aGaba_syn+0.5*K2(11)],p8,p17,p18,p19,p20,p21,Q.aI_syn,Q.aIgaba_syn,Q.ax,Q.G_syn,Q.Ggaba_syn);
K4=dt*Astrocyte(dt,[Q.ca+K3(1);Q.ap+K3(2);Q.aO1+K3(3);Q.aO2+K3(4);Q.aO3+K3(5); ...
                    Q.aR_syn+K3(6);Q.aE_syn+K3(7);Q.aRgaba_syn+K3(8);Q.aEgaba_syn+K3(9); ...
                    Q.aG_syn+K3(10);Q.aGaba_syn+K3(11)],p8,p17,p18,p19,p20,p21,Q.aI_syn,Q.aIgaba_syn,Q.ax,Q.G_syn,Q.Ggaba_syn);

Q.ca        = Q.ca        + (K1(1)+2*K2(1)+2*K3(1)+K4(1))/6;
Q.ap        = Q.ap        + (K1(2)+2*K2(2)+2*K3(2)+K4(2))/6;
Q.aO1       = Q.aO1       + (K1(3)+2*K2(3)+2*K3(3)+K4(3))/6;
Q.aO2       = Q.aO2       + (K1(4)+2*K2(4)+2*K3(4)+K4(4))/6;
Q.aO3       = Q.aO3       + (K1(5)+2*K2(5)+2*K3(5)+K4(5))/6;
Q.aR_syn    = Q.aR_syn    + (K1(6)+2*K2(6)+2*K3(6)+K4(6))/6;
Q.aE_syn    = Q.aE_syn    + (K1(7)+2*K2(7)+2*K3(7)+K4(7))/6;
Q.aRgaba_syn= Q.aRgaba_syn+ (K1(8)+2*K2(8)+2*K3(8)+K4(8))/6;
Q.aEgaba_syn= Q.aEgaba_syn+ (K1(9)+2*K2(9)+2*K3(9)+K4(9))/6;
Q.aG_syn    = Q.aG_syn    + (K1(10)+2*K2(10)+2*K3(10)+K4(10))/6;
Q.aGaba_syn = Q.aGaba_syn + (K1(11)+2*K2(11)+2*K3(11)+K4(11))/6;
Q.aI_syn    = 1 - Q.aR_syn - Q.aE_syn;
Q.aIgaba_syn= 1 - Q.aRgaba_syn - Q.aEgaba_syn;

Q.ax        = ax_next;

%% ===== Post-synaptic (NTS neuron) =====
g_ampa = p13(3); V_ampa=p13(4);
g_nmda = p14(3); V_nmda=p14(4);
g_gaba = p15(3); V_GABAa=p15(4);
g_gabab= p16(3); V_GABAb=p16(4);
Mg = p14(5);

I_postampa   = g_ampa*Q.m_ampa*(Q.V_post - V_ampa);
I_postgaba   = g_gaba*Q.m_gaba*(Q.V_post - V_GABAa);
I_postgabab  = g_gabab*Q.m_gaba*(Q.v - V_GABAb);
Mgv          = (1+exp((-0.062*Q.V_post*Mg)/3.57))^(-1);
I_postnmda   = g_nmda*Q.m_nmda*(Q.v - V_nmda)*Mgv;
I_postsum    = I_postampa + I_postnmda;

K1=dt*PostSynN(0,[Q.m_ampa;Q.m_nmda;Q.r_gaba;Q.r_gabab;Q.V_post], ...
               p12,p13,p14,p15,p16,Q.G_syn,Q.Ggaba_syn,Q.aGaba_syn,I_postsum,I_postgaba,I_postgabab);
K2=dt*PostSynN(0.5*dt,[Q.m_ampa+0.5*K1(1);Q.m_nmda+0.5*K1(2);Q.r_gaba+0.5*K1(3);Q.r_gabab+0.5*K1(4);Q.V_post+0.5*K1(5)], ...
               p12,p13,p14,p15,p16,Q.G_syn,Q.Ggaba_syn,Q.aGaba_syn,I_postsum,I_postgaba,I_postgabab);
K3=dt*PostSynN(0.5*dt,[Q.m_ampa+0.5*K2(1);Q.m_nmda+0.5*K2(2);Q.r_gaba+0.5*K2(3);Q.r_gabab+0.5*K2(4);Q.V_post+0.5*K2(5)], ...
               p12,p13,p14,p15,p16,Q.G_syn,Q.Ggaba_syn,Q.aGaba_syn,I_postsum,I_postgaba,I_postgabab);
K4=dt*PostSynN(dt,[Q.m_ampa+K3(1);Q.m_nmda+K3(2);Q.r_gaba+K3(3);Q.r_gabab+K3(4);Q.V_post+K3(5)], ...
               p12,p13,p14,p15,p16,Q.G_syn,Q.Ggaba_syn,Q.aGaba_syn,I_postsum,I_postgaba,I_postgabab);

Q.m_ampa = Q.m_ampa + (K1(1)+2*K2(1)+2*K3(1)+K4(1))/6;
Q.m_nmda = Q.m_nmda + (K1(2)+2*K2(2)+2*K3(2)+K4(2))/6;
Q.r_gaba = Q.r_gaba + (K1(3)+2*K2(3)+2*K3(3)+K4(3))/6;
Q.r_gabab= Q.r_gabab+ (K1(4)+2*K2(4)+2*K3(4)+K4(4))/6;
Q.V_post = Q.V_post + (K1(5)+2*K2(5)+2*K3(5)+K4(5))/6;

% Store currents (optional bookkeeping, not used elsewhere)
Q.I_postampa = I_postampa;
Q.I_postnmda = I_postnmda;
Q.I_postgaba = I_postgaba;
Q.I_postsum  = I_postsum;

% Advance absolute quadripartite time
Q.t_ms = Q.t_ms + dt;
end