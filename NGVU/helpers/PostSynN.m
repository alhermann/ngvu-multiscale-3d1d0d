% Quad-partite Synapse Model, Postsynaptic Neuron Dynamics
function [dydt] = PostSynN(~,y,param12,param13,param14,param15, ...
    param16,G_syn,Ggaba_syn,aGaba_syn,I_postsum,I_postgaba,I_postgabab)

% Constants and parameters
R_in=param12(1);
tau_mem = param12(2);
v_post=param12(6);
alpha_ampa=param13(1);
beta_ampa=param13(2);
alpha_nmda=param14(1);
beta_nmda=param14(2);
alpha_gaba=param15(1);
beta_gaba=param15(2);
alpha_gabab=param16(1);
beta_gabab=param16(2);

% Values set to equal input values
m_ampa = y(1);
m_nmda = y(2);
r_gaba = y(3);
r_gabab = y(4);
V_post = y(5);

% To simulate GABAb-Rs on the postsynaptic membrane, activate the I_postgabab
% current

% Postsynaptic Neuron response
dydt = [(alpha_ampa*G_syn*(1-m_ampa)-beta_ampa*m_ampa);                     % AMPA-R
    (alpha_nmda*G_syn*(1-m_nmda)-beta_nmda*m_nmda);                         % NMDA-R
    (alpha_gaba*Ggaba_syn*(1-r_gaba)-beta_gaba*r_gaba);                     % GABAa-R
    (alpha_gabab*(Ggaba_syn+aGaba_syn)*(1-r_gabab)-beta_gabab*r_gabab);     % GABAb-R
    (1/tau_mem)*(-(V_post-v_post)-R_in*(I_postsum+I_postgaba+I_postgabab))];

    end