function [lambda] = Switching(t,delta)
%Switching produces a smooth (infinitely differentiable) Switching function that
% switches on over time delta, remains constant, and then switches off over
% time delta.
%    INPUT:       t = time of interaction
%                 T = ramp-up time
%
%    OUTPUT: lambda = Switching function lambda(t)

dt                    = t(2)-t(1);
tf                    = t(length(t));
ndeltaON              = ceil(delta/dt);
ndeltaOFF             = length(t)-ndeltaON;
lambda                = ones(1,length(t));
lambda(1:ndeltaON)    = (1-tanh(cot(pi*t(1:ndeltaON)/delta)))/2;
lambda(ndeltaOFF:end) = (1-tanh(cot(pi*(tf-t(ndeltaOFF:length(t)))/delta)))/2;
end