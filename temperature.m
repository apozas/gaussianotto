%   Code used for creating Fig. 2 in arXiv:
%   Models the interaction between a harmonic oscillator machine and a ring
%   of harmonic oscillators initialized in a thermal state. After the
%   interaction the effective temperature of the machine is computed.
% 
%   authors:     Alejandro Pozas-Kerstjens, Karen V. Hovhanissyan,
%				 Eric G. Brown
%
%   requires:    -
%
%   last update: Aug, 2017
addpath('functions')
% -------------------------------------------------------------------------
% Choice of parameters
% -------------------------------------------------------------------------
% Machine frequency
Om = 2;

% Interaction strength with the bath
strength = 0.1;

% Bath and machine temperature
Tb = 4;
Tm = 0.5;

% Number of modes in the bath
N = floor(linspace(2,20,19));

% Time of interaction with the bath; must be at least 2*delta
tf = linspace(1,200,200);    

% Set of bath modes with which the machine interacts
interact = [1];

% -------------------------------------------------------------------------
% Computations
% -------------------------------------------------------------------------
% Initial detector state
sigDetI = eye(2)*(exp(Om/Tm)+1)/(exp(Om/Tm)-1);
    
Teff = zeros(length(tf),length(N)); % Matrix of temperatures

for n=1:length(N)
    for time=1:length(tf)
	    % Bath oscillators' frequencies
		freqs = Om*ones(1,N(n));
		% Baths' nearest-neighbour couplings
		alpha = 0.1*ones(1,N(n)-1);
		% Ramp-up time
		delta = 0.2*tf(time);
		% Bath free Hamiltonian.
        Ffree = FreeRing(N(n),freqs,alpha);
        % Initialize global (machine+bath) state
        sigI = blkdiag(sigDetI,Initialize(N(n),Tb,Ffree));
        %Initialize time steps for numerical integrations
        dt      = 0.01/Om;
		steps   = floor(tf(time)/dt);
        dt      = tf(time)/steps;    % Recompute to account for rounding
        t       = linspace(0,tf(time),steps);
        lambda  = Switching(t,delta);
        S       = MakeStimeIndep(N(n),Om,strength,interact,t,delta,lambda,Ffree);
        sigF    = S*sigI*S';
        sigDetF = sigF(1:2,1:2);
		eigen   = sort(eigs(sigDetF));
        nu      = sqrt(eigen(1)*eigen(2));
        
        Teff(time,n) = Om/(log((nu+1)/(nu-1)));
    end
end

% Temperature vs. size and interaction time (Figure 2)
figure(1)
imagesc(N,tf,Teff)
set(gca,'Ydir','normal')
colorbar
title('T_{eff}')
xlabel('N')
ylabel('\tau')