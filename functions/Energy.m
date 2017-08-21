function [E] = Energy(sig,Ffree,Om)
%Energy computes the free energy of the bath-machine system
%    INPUT:   sig = machine-bath(s) covariance matrix
%           Ffree = Free machine-bath(s) Hamiltonian matrix
%              Om = Machine's frequency
%
%    OUTPUT:    E = system's energy

Ffree = blkdiag(Om*eye(2)/2,Ffree);
E     = trace(Ffree*sig)/2;
end