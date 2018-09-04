%   Code used for creating Figs. 5 and 6 in New J. Phys. 20, 043034 (2018)
%   Runs an Otto cycle between two baths of harmonic oscillators, computing
%   the mutual information between the machine oscillator and each of the
%   oscillators in the baths, the mutual information between the machine
%   oscillator and the complete bath, and between every pair of oscillators
%   in the baths, for a number of interaction cycles.
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
% Machine frequencies
OmC       = 1;
OmH       = 2;

% Machine coupling strengths with the baths
gammaC    = 0.1;    
gammaH    = 0.1;

% Baths temperatures
TC        = 0.5;
TH        = 4;

% Number of modes in the baths
NC        = 30;
NH        = NC;

% Baths modes' frequencies
freqsC    = OmC*ones(1,NC);
freqsH    = OmH*ones(1,NH);

% Baths inner couplings (nearest-neighbors)
alphaC    = 0.1*ones(1,NC-1);
alphaH    = 0.1*ones(1,NH-1);

% Number of cycles we run for
Ncycle    = 10;

% Time of interaction with one bath; must be at least 2*delta
tf        = 100;    

% Ramp-up time; want it to be large compared to inverse detector 
% frequency so as to remain approximately adiabatic
delta     = 0.1*tf;

% Time steps for numerical integrations; must be very small
% compared to inverse detector frequency
dt        = 0.01;

% Set of bath modes with which the machine interacts
interactC = [floor(NC/2)];
interactH = NC + interactC;

% -------------------------------------------------------------------------
% Computations
% -------------------------------------------------------------------------
% Baths free Hamiltonians
FfreeC   = FreeRing(NC,freqsC,alphaC);
FfreeH   = FreeRing(NH,freqsH,alphaH);

% Initial detector state
% We start by interacting with the hot bath, assuming that the WM
% is in a state thermal with the cold bath
sigDetI   = eye(2)*(exp(OmC/TC)+1)/(exp(OmC/TC)-1);

% Initialize global state of the system (detector plus two baths)
sigI      = blkdiag(sigDetI,Initialize(NC,TC,FfreeC),Initialize(NH,TH,FfreeH));

% Switching function
steps     = floor(tf/dt);
dt        = tf/steps;    % Recompute dt to account for rounding
t         = linspace(0,tf,steps);
lambda    = Switching(t,delta);

% Compute the symplectic evolution matrix for both bath interactions.
SC = MakeS(NC+NH,OmC,gammaC,interactC,t,delta,lambda,blkdiag(FfreeC,FfreeH));
SH = MakeS(NC+NH,OmH,gammaH,interactH,t,delta,lambda,blkdiag(FfreeC,FfreeH));

% Computation of mutual informations.
% We look at quantities every 100 time steps.
ntH           = floor(length(t)/100);
ntC           = ntH;
MutInf        = zeros(NC+NH,ntH+ntC,Ncycle);
MutInfTotal   = zeros(Ncycle,ntH+ntC);
MutInfLattice = zeros(NC+NH,NC+NH,(ntH+ntC)*Ncycle);

for cyc=1:Ncycle
    % Interaction with the hot bath
    for i=1:ntH
        % State evolution
        sigH      = SH(:,:,100*i-99)*sigI*SH(:,:,100*i-99)';
		sigDetH   = sigH(1:2,1:2);
		sigBathsH = sigH(3:end,3:end);
        % Total mutual information between WM and baths
        MutInfTotal(cyc,i) = Entropy(sigDetH)+Entropy(sigBathsH) ...
                            -Entropy(sigH);
        % Mutual information between WM and individual bath modes 
        for j=1:NC+NH
		    sigMode = sigH(2*j+1:2*j+2,2*j+1:2*j+2);
            sigComb = [sigDetH,sigH(1:2,2*j+1:2*j+2);sigH(2*j+1:2*j+2,1:2),sigMode];
            MutInf(j,i,cyc) = Entropy(sigDetH)+Entropy(sigMode) ...
                             -Entropy(sigComb);
            for k=1:NC+NH
                % Mutual information between each pair of bath modes
                if (j ~= k)
				    sigMode2 = sigH(2*k+1:2*k+2,2*k+1:2*k+2);
                    sigComb  = [sigMode,sigH(2*j+1:2*j+2,2*k+1:2*k+2);sigH(2*k+1:2*k+2,2*j+1:2*j+2),sigMode2];
                    MutInfLattice(j,k,(ntH+ntC)*(cyc-1)+i) = Entropy(sigMode)+Entropy(sigMode2)...
                                                            -Entropy(sigComb);
                end
            end
        end
    end
	% Set state of the system after the interaction with the hot bath
    sigI = SH(:,:,end)*sigI*SH(:,:,end)';
    % Interaction with the cold bath
    for i=1:ntC
        % State evolution
        sigC      = SC(:,:,100*i-99)*sigI*SC(:,:,100*i-99)';
		sigDetC   = sigC(1:2,1:2);
		sigBathsC = sigC(3:end,3:end);
        % Total mutual information between WM and baths
        MutInfTotal(cyc,i+ntH) = Entropy(sigDetC)+Entropy(sigBathsC) ...
                                -Entropy(sigC);
        % Mutual information between WM and individual bath modes 
        for j=1:NC+NH
		    sigMode = sigC(2*j+1:2*j+2,2*j+1:2*j+2);
            sigComb = [sigDetC(1:2,1:2),sigC(1:2,2*j+1:2*j+2);sigC(2*j+1:2*j+2,1:2),sigMode];
            MutInf(j,i+ntH,cyc) = Entropy(sigDetC)+Entropy(sigMode) ...
                                 -Entropy(sigComb);
            for k=1:NC+NH
                % Mutual information between each pair of bath modes
                if (j ~= k)
				    sigMode2 = sigC(2*k+1:2*k+2,2*k+1:2*k+2);
                    sigComb  = [sigMode,sigC(2*j+1:2*j+2,2*k+1:2*k+2);sigC(2*k+1:2*k+2,2*j+1:2*j+2),sigMode2];
                    MutInfLattice(j,k,(ntH+ntC)*(cyc-1)+ntH+i) = Entropy(sigMode)+Entropy(sigMode2) ...
                                                                -Entropy(sigComb);
                end
            end
        end
    end
	% Set state of the system after the interaction with the cold bath
    sigI = SC(:,:,end)*sigI*SC(:,:,end)';
end

% -------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------
% Total mutual information
figure(1)
plot(linspace(1,Ncycle*2*t(end),Ncycle*(ntC+ntH)),reshape(MutInfTotal',Ncycle*(ntC+ntH),1))
y1=get(gca,'ylim');
hold on
for i=1:Ncycle-1
    plot([2*t(end)*i 2*t(end)*i],y1,'Color','k')
    plot([t(end)*(2*i-1) t(end)*(2*i-1)],y1,'Color','k', 'LineStyle', '--')
end
plot([t(end)*(2*Ncycle-1) t(end)*(2*Ncycle-1)],y1,'Color','k', 'LineStyle', '--')
hold off
title('Machine-Baths total mutual information')

% WM-bath modes MI (Figure 5)
MutInfLong = zeros(NC+NH,(ntC+ntH)*Ncycle);
for i=1:Ncycle
    MutInfLong(:,(ntC+ntH)*(i-1)+1:(ntC+ntH)*i) = MutInf(:,:,i);
end
figure(2)
imagesc(linspace(1,Ncycle*2*t(end),Ncycle*(ntC+ntH)),linspace(1,NC+NH,NC+NH),MutInfLong)
colormap('hot')
colorbar
y1=get(gca,'ylim');
hold on
for i=1:Ncycle-1
    plot([2*t(end)*i 2*t(end)*i],y1,'Color','w')
    plot([t(end)*(2*i-1) t(end)*(2*i-1)],y1,'Color','w', 'LineStyle', '--')
end
plot([t(end)*(2*Ncycle-1) t(end)*(2*Ncycle-1)],y1,'Color','w', 'LineStyle', '--')
hold off
title('Machine-Baths mutual information')

% Animation of mutual information between pairs of oscillators (Figure 6)
figure(3)
for i=1:length(MutInfLattice(1,1,:))
    figure(3)
    imagesc(linspace(1,NC+NH,NC+NH),linspace(1,NC+NH,NC+NH),MutInfLattice(:,:,i))
    colormap('hot')
    colorbar
    title('Bath-Bath mutual information')
    ylabel('    Hot bath                                   Cold bath')
    xlabel('Cold bath                                       Hot bath')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    pause(0.01)
end