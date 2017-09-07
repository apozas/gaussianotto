function [Fint] = MakeInt(N,interact)
%MakeInt constructs the interaction part of the Hamiltonian matrix.
%    Called by MakeStimeIndep
%    INPUT:         N = number of modes in the bath
%            interact = set of modes the WM interacts with
%
%    OUTPUT:     Fint = Interaction Hamiltonian matrix

% Make interaction submatrix
X               = zeros(1,2*N); % Row of interactions with the machines' q
X(2*interact-1) = 1/2;
X               = vertcat(X,zeros(1,2*N)); % Add row of interactions with the machines' p

% X is in the off-diagonal blocks of Fint
Fint = horzcat(zeros(2),X);
Fint = vertcat(Fint,horzcat(X',zeros(2*N)));
end