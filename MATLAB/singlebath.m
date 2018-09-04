%   Code used for creating Figs. 1 and 7 in New J. Phys. 20, 043034 (2018)
%   Models the interaction between a harmonic oscillator WM and a ring
%   of harmonic oscillators initialized in a thermal state. During the
%   interaction the temperature of the WM, its athermality, and three
%   mutual informations (WM-bath, WM-interacting oscillator and
%   interacting oscillator-rest of the bath) are computed.
%
%   authors:     Alejandro Pozas-Kerstjens, Karen V. Hovhanissyan,
%                Eric G. Brown
%
%   requires:    -
%
%   last update: Sep, 2017
addpath('functions')
% -------------------------------------------------------------------------
% Choice of parameters
% -------------------------------------------------------------------------
% WM frequency
Om       = 2;

% Interaction strength with the bath
gamma    = 0.1;

% Bath and machine temperature
Tb       = 4;
Tm       = 0.5;

% Number of modes in the bath
N        = 30;

% Bath oscillators' frequencies
freqs    = Om*ones(1,N);

% Baths' nearest-neighbour couplings
alpha    = 0.1*ones(1,N-1);

% Time of interaction with the bath; must be at least 2*delta
tf       = 100;    

% Ramp-up time; want it to be large compared to inverse WM's 
% frequency so as to remain approximately adiabatic
delta    = 0.1*tf;

% Time steps for numerical integrations; must be very small
% compared to inverse WM's frequency
dt       = 0.01;

% Set of bath modes with which the WM interacts
interact = [1];

% -------------------------------------------------------------------------
% Computations
% -------------------------------------------------------------------------
% Bath free Hamiltonian.
Ffree   = FreeRing(N,freqs,alpha);

% Initial detector state
sigDetI = eye(2)*(exp(Om/Tm)+1)/(exp(Om/Tm)-1);

% Initialize global state of the system (WM plus bath)
sigI    = blkdiag(sigDetI,Initialize(N,Tb,Ffree));

% Switching function
steps   = floor(tf/dt);
dt      = tf/steps;    % Recompute dt to account for rounding
t       = linspace(0,tf,steps);
lambda  = Switching(t,delta);

% Compute the symplectic evolution for both bath interactions.
S       = MakeS(N,Om,gamma,interact,t,delta,lambda,Ffree);

% Computation of evolution of the WM's temperature and its athermality
Teff    = zeros(1,length(t));
Teff(1) = Tm;
ath     = zeros(1,length(t));
for i=2:length(t)
    sigF    = S(:,:,i)*sigI*S(:,:,i)';
    sigDetF = sigF(1:2,1:2);
    eigen   = sort(eigs(sigDetF));
    nu      = sqrt(eigen(1)*eigen(2));
    Teff(i) = Om/(log((nu+1)/(nu-1)));
    thermal = eye(2)*(exp(Om/Teff(i))+1)/(exp(Om/Teff(i))-1);
    Delta   = 4*det(sigDetF+thermal);
    Lambda  = (4*det(thermal)-1)*(4*det(sigDetF)-1);
    ath(i)  = 1-2/(sqrt(Delta+Lambda)-sqrt(Lambda));
end

% Computation of mutual informations.
% We look at quantities every 100 time steps.
nt            = floor(length(t)/100);
MutInfdetbath = zeros(nt,1); % MI between WM and bath
MutInfdetosc  = zeros(nt,1); % MI between WM and interacting mode
MutInfoscrest = zeros(nt,1); % MI between interacting mode and the rest

for i=1:nt
    % Compute evolved state
    sig = S(:,:,100*i-99)*sigI*S(:,:,100*i-99)';
    
    % Total mutual information between WM and bath
    MutInfdetbath(i) = Entropy(sig(1:2,1:2))+Entropy(sig(3:end,3:end)) ...
                      -Entropy(sig);
    % Total mutual information between WM and interacting oscillator
    sigComb          = sig(1:4,1:4);
    MutInfdetosc(i)  = Entropy(sig(1:2,1:2))+Entropy(sig(3:4,3:4)) ...
                      -Entropy(sigComb);
    % Total mutual information between interacting oscillator and the rest
    sigComb          = sig(3:end,3:end);
    MutInfoscrest(i) = Entropy(sig(3:4,3:4))+Entropy(sig(5:end,5:end)) ...
                      -Entropy(sigComb);
end

% -------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------
% Various correlations and temperature (Figure 1)
figure(1)
plot(linspace(1,t(end),nt),MutInfdetbath)
hold on
plot(linspace(1,t(end),nt),MutInfdetosc,'k')
plot(linspace(1,t(end),nt),MutInfoscrest,'r')
plot(t,Teff,'g-')
hold off
legend('Machine-Bath','Machine-Osc','Osc-Rest','T_{eff}','Location','east')
xlabel('t')

% Athermality (Figure 7)
figure(2)
plot(t,ath)
ylabel('Athermality')
xlabel('t')