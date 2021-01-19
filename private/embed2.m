function Y2 = embed2(varargin)
% EMBED2 Extend a phase space vector by a new component.
%    Y2 = EMBED2(Y, X, TAU) extends the N-by-M matrix Y by a new
%    column using the time shifted values of time series X.
%    The time shift is specified by value TAU. The final matrix 
%    Y2 has size (N-tau)-by-(M+1).
 
% Copyright (c) 2020
% K. Hauke Kraemer, 
% Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software and runs under MIT licence.

%%
Y = varargin{1};
x = varargin{2};
tau = varargin{3};

if size(Y,1) < size(Y,2)
    Y = Y';
end

if size(x,1) < size(x,2)
    x = x';
end

N = size(Y,1);
MM = length(x);
MMM = MM - tau;
M = min([N, MMM]);
Y2 = [Y(1:M, :) x(1+tau:tau+M)];
