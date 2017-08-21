function [S] = MakeS(N,Om,strength,interact,t,delta,lambda,Ffree)
%MakeS constructs the symplectic evolution matrix of an interaction.
%    The output is given in an array where the matrix in position (i) is the
%    time evolution matrix from the beginning until t(i)
%    INPUT:         N = number of modes in the bath
%                  Om = modes frequency
%            strength = bath-machine coupling strength
%            interact = set of modes the machine interacts with
%                   t = array of time points used for numerical integration
%               delta = ramp-up time
%              lambda = array of values of the Switching at the steps in t
%               Ffree = Bath free Hamiltonian matrix
%
%    OUTPUT:        S = Array of time evolution matrices

Form  = kron(eye(1+N),[0,1;-1,0]);    % Symplectic form in p-q basis

Ffree = blkdiag(Om*eye(2)/2,Ffree);   % Free Hamiltonian matrix
Fint  = MakeInt(N,interact);          % Interaction Hamiltonian matrix

dt        = t(2)-t(1);
ndeltaON  = ceil(delta/dt);
ndeltaOFF = length(t)-ndeltaON;

Sdt = zeros(2*(N+1),2*(N+1),ndeltaON+1); % Minimal array of time evolution matrices

% Note that the Hamiltonians involved are symmetric, and thus we compute
% 2*F instead of F+F^T
for i=1:ndeltaON+1
    Sdt(:,:,i) = expm(2*Form*(Ffree+strength*lambda(i)*Fint)*dt);
end

S        = zeros(2*(N+1),2*(N+1),length(t));    % Time evolution array
S(:,:,1) = Sdt(:,:,1);
% Ramp-up
for i=2:ndeltaON
    S(:,:,i) = Sdt(:,:,i)*S(:,:,i-1);
end
% Plateau
Splat=expm(2*Form*(Ffree+strength*Fint)*dt); % Time evolution in the plateau for time dt
for i=ndeltaON+1:ndeltaOFF
    S(:,:,i) = Splat*S(:,:,i-1);
end
% Ramp-down
for i=ndeltaOFF+1:length(t)
    S(:,:,i) = Sdt(:,:,length(t)-(i-1))*S(:,:,i-1);
end
end