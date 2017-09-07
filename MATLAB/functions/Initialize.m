function [sigBathI] = Initialize(N,T,Ffree)
%Initialize computes covariance matrix of bath at temperature T, given the local
%    free Hamiltonian matrix Ffree
%
%    NOTE: This algorithm only works assuming that the free Hamiltonian
%    contains no q-p interaction terms. q-q and p-p terms are valid.
%
%    INPUT:         N = number of oscillators in the bath
%                   T = Bath temperature
% 		        Ffree = Free bath Hamiltonian matrix
%
%    OUTPUT: sigBathI = Bath covariance matrix

% First transform Ffree to q-p Ordering
Z     = Ordering(N);
Ffree = Z'*Ffree*Z;

% Computation of transformation matrix to normal-mode basis
Fq        = Ffree(1:N,1:N);
Fp        = Ffree(N+1:2*N,N+1:2*N);
A         = Fq^(1/2)*Fp^(1/2);
[O1,d,O2] = svd(A);
O1        = O1';
O2        = O2';                                              
O         = blkdiag(O1,O2);
S         = blkdiag(d,d)^(1/2)*O*Ffree^(-1/2);
S         = S';                   % Transformation matrix

% Construct the state in the normal basis
normals   = 2*svd(A);         % Normal frequencies
nu        = (exp(normals/T)+1)./(exp(normals/T)-1);    % Thermal symplectic eigenvalues
sigNormal = blkdiag(diag(nu),diag(nu));

% Transformation to local basis
sigLocal = S*sigNormal*S';

% Transform back to mode-mode Ordering
sigBathI = Z*sigLocal*Z';
end