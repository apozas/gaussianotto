function [E] = Energy(sig,Ffree,Om)
%Energy computes the free energy of the bath(s)-WM system
%    INPUT:   sig = WM-bath(s) covariance matrix
%           Ffree = Bath(s) free Hamiltonian matrix
%              Om = WM's frequency
%
%    OUTPUT:    E = System's energy

Ffree = blkdiag(Om*eye(2)/2,Ffree);
E     = trace(Ffree*sig)/2;
end