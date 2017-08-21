function [H] = Entropy(sig)
%Entropy computes Entropy of a given covariance matrix
%    INPUT: sig = covariance matrix
%
%    OUTPUT:  H = Entropy
siz    = size(sig);
N      = siz(1)/2;
Form   = kron(eye(N),[0,1;-1,0]);    % Symplectic form in mode-mode basis
nufull = sort(abs(eig(j*Form*sig))); % Eigenvalues in +-nu pairs
nu     = zeros(N,1);
for m=1:N
    nu(m) = nufull(2*m-1);
end

H = sum((nu+1).*log((nu+1)/2)/2-(nu-1).*log((nu-1)/2)/2);
end