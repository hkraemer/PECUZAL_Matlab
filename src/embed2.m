function Y2 = embed2(varargin)
% EMBED2 takes a matrix 'Y' containing all phase space vectors, a univariate
% time series 'x' and a tau value 'tau' as input. embed2 then expands the 
% input phase space vectors by an additional component consisting of the 
% tau-shifted values of the input time series x.
% 
%               Y2 = embed2(Y,x,tau)
% 
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

if size(Y,1)<size(Y,2)
    Y = Y';
end

if size(x,1)<size(x,2)
    x = x';
end

N = size(Y,1); 

timespan_diff = tau;
M = N - timespan_diff;

Y2 = zeros(M,size(Y,2)+1);
Y2(:,1:size(Y,2)) = Y(1:M,:);

Y2(:,size(Y,2)+1) = x(1+tau:N);

end

