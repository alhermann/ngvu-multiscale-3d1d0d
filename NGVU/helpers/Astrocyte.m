% Quad-partite Synapse Model, Actrocyte Dynamics
function [dydt] = Astrocyte(~,y,param8,param17,param18,param19, ...
            param20,param21,aI_syn,aIgaba_syn,ax,G_syn,Ggaba_syn)

% Constants and parameters
np = param8(6);
ac0=param17(2);
ac1=param17(3);
av1=param17(4);
av2=param17(5);
av3=param17(6);
ak3=param17(7);
ad1=param17(8);
ad5=param17(11);
av_plcd=param18(1);
ak_plcd=param18(2);
aK_plcd=param18(3);
ar5p=param18(4);
av_plcb=param18(5);
avgaba_plcb=param18(6);
aK_R=param18(7);
aK_R_gaba=param18(8);
aK_P=param18(9);
aK_P_gaba=param18(10);
aK_pi=param18(11);
aK_pi_gaba=param18(12);
av_3K=param18(13);
aK_D=param18(14);
aK3=param18(15);
nva=param19(1);
nva_gaba=param19(2);
gva=param19(3);
gva_gaba=param19(4);
ak1=param19(5);
ak2=param19(6);
akk3=param19(7);
ak_1=param19(8);
ak_2=param19(9);
ak_3=param19(10);
adegG=param20(1);
atau_rec=param20(2);
atau_inact=param20(3);
adeg_gaba=param21(1);
atau_rec_gaba=param21(2);
atau_inact_gaba=param21(3);

% Values set to equal input values
ca = y(1);
ap = y(2);
aO1 = y(3);
aO2 = y(4);
aO3 = y(5);
aR_syn= y(6);
aE_syn = y(7);
aRgaba_syn = y(8);
aEgaba_syn=y(9);
aG_syn=y(10);
aGaba_syn=y(11);

% Time-constants for three binding sites on the SLMV
atau1=(ak1*ca+ak_1);                                                % Closure of S1 (1/ms)
atau2=(ak2*ca+ak_2);                                                % Closure of S2 (1/ms)
atau3=(akk3*ca+ak_3);                                               % Closure of S3 (1/ms)

% IP3 production & degradation by astrocyte
aplcb_ca=1+(aK_P/aK_R)*(ca/(ca+aK_pi));                             % Ca++-dependent inhibition of 'ap_glu'
ap_glu=av_plcb*(G_syn^(np-.4))/((G_syn^(np-.4))+ ...
    (aK_R*aplcb_ca)^(np-.4));                                       % Glu-dependent IP3 production
aplcb_ca_gaba=1+(aK_P_gaba/aK_R_gaba)*(ca/(ca+aK_pi_gaba));         % Ca++-dependent inhibition of 'ap_gaba'
ap_gaba=avgaba_plcb*(Ggaba_syn^(np-.4))/((Ggaba_syn^(np-.4))+ ...
    (aK_R_gaba*aplcb_ca_gaba)^(np-.4));                             % GABA-dependent IP3 production
ap_plcd=av_plcd*(1/(1+(ap/ak_plcd)))*...
    (ca^2)/((ca^2)+aK_plcd^2);                                      % Neurotransmitter-independent IP3 production
ap_mapk=av_3K*((ca^4)/((ca^4)+ ...
    aK_D^4))*(ap/(ap+aK3));                                         % IP3 degradation by IP3-3K
ap_deg=ar5p*(ap);                                                   % IP3 degradation by IP-5P

% IP3-R dynamics
aminf=ap/(ap+ad1);
aninf=ca/(ca+ad5);                                                  % Steady-state functions for IP3 & Ca2+
auer=(ac0-ca)/ac1;
ajchan=ac1*av1*(aminf^3)*(aninf^3)*((abs(ax))^3)*(ca-auer);         % Ca++ flux through IP3-R
ajpump=av3*(ca^2)/(ak3^2+ca^2);                                     % SERCA pump
ajleak=ac1*av2*(ca-auer);                                           % Ca++ leak from the cytosol to the ER

% Astrocyte Dynamics
dydt = [(-ajchan-ajpump-ajleak);                                                  % [Ca++]astr turnover
    (ap_glu+ap_gaba+ap_plcd-ap_mapk-ap_deg);                                      % The total amount of [IP3] produced
    (ak1*ca-(aO1*atau1));                                                         % Site 1:SLMV release Ca++ bound
    (ak2*ca-(aO2*atau2));                                                         % Site 2:SLMV release Ca++ bound
    (akk3*ca-(aO3*atau3));                                                        % Site 3:SLMV release Ca++ bound
    (((aI_syn)/atau_rec)-((heaviside(ca-196.69)*aO1*aO2*aO3)*aR_syn));            % Fraction of Glu containing releasable SLMVs
    (((heaviside(ca-196.69)*aO1*aO2*aO3)*aR_syn)-(aE_syn/atau_inact));            % Fraction of Glu containing effective SLMVs
    (((aIgaba_syn)/atau_rec_gaba)-((heaviside(ca-196.69)*aO1*aO2*aO3)*aRgaba_syn));    % Fraction of GABA containing releasable SLMVs
    (((heaviside(ca-196.69)*aO1*aO2*aO3)*aRgaba_syn)-(aEgaba_syn/atau_inact_gaba));    % Fraction of GABA containing effective SLMVs
    (nva*gva*aE_syn-adegG*aG_syn);                                              % [Glu] released from the astrocyte (nM)
    (nva_gaba*gva_gaba*aEgaba_syn-adeg_gaba*aGaba_syn)];                        % [GABA] released from the astrocyte (nM);       
end