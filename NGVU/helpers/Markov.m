% Markov.m predicts the present state of the Ca++ sensor depending upon its 
% previous state. It generates a URN 'u' and returns index 'i' for which u  
% is greater than p(i), where p is the transition probability matrix with
% the sum of each elements of the column = '1'. The method of the Markov
% process simulation is described in the book:
% 'Computational Cell Biology' by Christopher P. Fall et al (2002) 
% page number: 291-292.
% The value of i=1,2,3,...,7. If i=1 it denotes sensor has no Ca++ is 
% bound; if i=7 then the sensor has 5 sCa++ bound and is in the isomerized
% state (Equation (7) in the manuscript by Tewari, S. and Majumdar, K. 
% (2012). "A Mathematical Model of Tripartite Synapse: 
% Astrocyte Induced Synaptic 
% Plasticity," Journal of Biological Physics, 38(3): pp. 465–496.

function [state] = Markov(p)
u=rand;                         % URN
i=1;                            % Initial value of i
s=p(1);                         % Probability to remain in the current state

while ((u>s) && (i<length(p)))
  i=i+1; 
  s=s+p(i);
end

state=i;                        % i indicates the present state. 

