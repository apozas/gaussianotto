function [S] = MakeStimeIndep(N,Om,strength,interact,t,delta,lambda,Ffree)
%MakeStimeIndep constructs the symplectic evolution matrix of a bath-WM system
%    during an interaction.
%    INPUT:         N = number of modes in the bath
%                  Om = WM's frequency
%            strength = bath-WM coupling strength
%            interact = set of modes the WM interacts with
%                   t = array of time points used for numerical integration
%               delta = ramp-up time
%              lambda = array of values of Switching at the steps in t
%               Ffree = Bath free Hamiltonian matrix
%
%    OUTPUT:        S = Time evolution matrix

Form  = kron(eye(1+N),[0,1;-1,0]);    % Symplectic form in p-q basis
Ffree = blkdiag(Om*eye(2)/2,Ffree);   % WM-bath free Hamiltonian matrix
Fint  = MakeInt(N,interact);

dt        = t(2)-t(1);
ndeltaON  = ceil(delta/dt);
ndeltaOFF = length(t)-ndeltaON;

% Computation of S
% Evolution during the plateau. Since F is constant and lambda as well,
% integration is trivial.
% Note that the Hamiltonians involved are symmetric, and thus we compute
% 2*F instead of F+F^T
S = expm(2*Form*(Ffree+strength*Fint)*dt*length(t(ndeltaON+1:ndeltaOFF)));

% Evolution during ramp-up and down. It is symmetric going outwards from the plateau
for i=ndeltaOFF+1:length(t)
	tempS = expm(2*Form*(Ffree+strength*lambda(i)*Fint)*dt);
    S     = tempS*S*tempS;
end
end