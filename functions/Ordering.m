function [Z] = Ordering(N)
%Ordering produces the 2N*2N transformation between the q-p Ordering and the
%    mode-mode Ordering
%    INPUT:      N = number of modes in the system
%
%    OUTPUT: Ffree = free Hamiltonian matrix for the chain

ZL = zeros(2*N,N);
ZR = ZL;

for j=1:N
    ZL(2*j-1,j) = 1;
    ZR(2*j,j)   = 1;
end

Z = horzcat(ZL,ZR);
end