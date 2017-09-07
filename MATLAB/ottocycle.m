%   Code used for creating Fig. 4 in arXiv:1708.06363
%   Runs an Otto cycle between two baths of harmonic oscillators, computing
%   the work output after every interaction, the heat loss at every cycle,
%   and the efficiency of the machine after every cycle.
%   Variables with a C or H at the end refer to the parameters used for the
%   hot bath (H) and the cold bath (C).
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
% WM frequencies
OmC       = 1;
OmH       = 2;

% WM coupling strengths with the baths
gammaC    = 0.1;    
gammaH    = 0.1;

% Baths' temperatures
TC        = 0.5;
TH        = 4;

% Number of modes in the baths
NC        = 300;
NH        = NC;

% Baths modes' frequencies
freqsC    = OmC*ones(1,NC);
freqsH    = OmH*ones(1,NH);

% Baths inner couplings (nearest-neighbors)
alphaC    = 0.1*ones(1,NC-1);
alphaH    = 0.1*ones(1,NH-1);

% Number of cycles we run for
Ncycle    = 200;

% Time of interaction with one bath; must be at least 2*delta
tf        = 100;    

% Ramp-up time; want it to be large compared to inverse WM's 
% frequency so as to remain approximately adiabatic
delta     = 0.1*tf;

% Time steps for numerical integrations; must be very small
% compared to inverse WM's frequency
dt        = 0.01;

% Set of bath modes with which the WM interacts
interactC = [1];
interactH = NC + interactC;

% -------------------------------------------------------------------------
% Computations
% -------------------------------------------------------------------------
% Baths free Hamiltonians
FfreeC   = FreeRing(NC,freqsC,alphaC);
FfreeH   = FreeRing(NH,freqsH,alphaH);

% Initial detector state
% We start by interacting with the hot bath, assuming that the WM is in a
% state thermal with the cold bath
sigDetI  = eye(2)*(exp(OmC/TC)+1)/(exp(OmC/TC)-1);

% Initialize global state of the system (detector plus two baths)
sigI     = blkdiag(sigDetI,Initialize(NC,TC,FfreeC),Initialize(NH,TH,FfreeH));

% Switching function
steps    = floor(tf/dt);
dt       = tf/steps;    % Recompute dt to account for rounding
t        = linspace(0,tf,steps);
lambda   = Switching(t,delta);

% Compute the symplectic evolution matrix for both bath interactions.
SC = MakeStimeIndep(NC+NH,OmC,gammaC,interactC,t,delta,lambda,blkdiag(FfreeC,FfreeH));
SH = MakeStimeIndep(NC+NH,OmH,gammaH,interactH,t,delta,lambda,blkdiag(FfreeC,FfreeH));

% Begin cyclic operation
workAdiabatic   = zeros(1,2*Ncycle);  % Work output during adiabats
workInteraction = zeros(1,2*Ncycle);  % Work output during isochores
workcycle       = zeros(1,Ncycle);    % Work output per cycle
workTotal       = zeros(1,Ncycle);    % Cumulative total work output
Q               = zeros(1,Ncycle);    % Energy loss in the hot bath
efficiency      = zeros(1,Ncycle);    % Cycle efficiency
sig             = sigI;

for i=1:Ncycle
    sigPre = sig;        % Store pre-interaction state
	% Hot bath interaction
    sig                    = SH*sig*SH';        
    sigDet                 = sig(1:2,1:2);
    workInteraction(2*i-1) = Energy(sigPre,blkdiag(FfreeC,FfreeH),OmH) ...
                            -Energy(sig,blkdiag(FfreeC,FfreeH),OmH);
    workAdiabatic(2*i-1)   = (OmH-OmC)*trace(sigDet)/4;
	sigHPre                = sigPre(2*(NC+1)+1:end,2*(NC+1)+1:end);
	sigH                   = sig(2*(NC+1)+1:end,2*(NC+1)+1:end);
    Q(i)                   = trace(FfreeH*sigHPre)/2-trace(FfreeH*sigH)/2;
    
    sigPre = sig;
	% Cold bath interaction
    sig                    = SC*sig*SC';
    sigDet                 = sig(1:2,1:2);
    workInteraction(2*i)   = Energy(sigPre,blkdiag(FfreeC,FfreeH),OmC) ...
                            -Energy(sig,blkdiag(FfreeC,FfreeH),OmC);
    workAdiabatic(2*i)     = (OmC-OmH)*trace(sigDet)/4;
    workTotal(i)           = sum(workAdiabatic)+sum(workInteraction);
end

for i=1:Ncycle
    if i==1
        workcycle(i) = workTotal(i);
    else
        workcycle(i) = workTotal(i) - workTotal(i-1);
    end
	efficiency(i)    = workcycle(i) / Q(i);
end

% -------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------
% Work output in each adiabat (Figures 4a, 4b)
figure(1)
bar(workAdiabatic)
xlabel('Interaction')
ylabel('Work per adiabat')

% Cumulative work output
figure(2)
plot(workTotal)
xlabel('Cycle')
ylabel('Cumulative work')

% Work output in each cycle (Figures 4a, 4b, 4c)
figure(3)
plot(workcycle)
xlabel('Cycle')
ylabel('Work per cycle')

% Energy loss of the hot bath (Figure 4c)
figure(4)
plot(Q)
xlabel('Cycle')
ylabel('Q')

% Machine's efficiency (Figure 4c)
figure(5)
plot(efficiency)
xlabel('Cycle')
ylabel('Efficiency')