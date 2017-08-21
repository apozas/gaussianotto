function [Ffree] = FreeRing(N,freqs,alpha)
%FreeRing constructs the free Hamiltonian matrix for an oscillator chain of the
%    form H = 1/2 sum_n (freq(n)*(q_n^2+p_n^2)+alpha(n)*q_n*q_(n+1), in the
%    mode-mode basis
%    INPUT:      N = number of oscillators in the chain
%            freqs = frequency of each oscillator
% 		     alpha = coupling strength between each pair of neighbor oscillators
%
%    OUTPUT: Ffree = free Hamiltonian matrix for the chain

Fq      = 0.5*(diag(freqs)+diag(alpha,1)+diag(alpha,-1));
Fq(1,N) = 0.5*alpha(1);
Fq(N,1) = 0.5*alpha(1);
Fp      = 0.5*diag(freqs);
F       = blkdiag(Fq,Fp);

% For consistency, transform to mode-mode Ordering
Z       = Ordering(N);
Ffree   = Z*F*Z';
end