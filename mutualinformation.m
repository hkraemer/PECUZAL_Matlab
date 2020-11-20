function Z = mutualinformation(varargin)
% MUTUALINFORMATION computes the automutualinformation values of timeseries.
%    I = MUTUALINFORMATION(X) returns the auto mutual information I of X,
%    where X can be a multi-column vector.
%
%    I = MUTUALINFORMATION(X, MAXTAU) returns the auto mutual information for
%    the lags in the range [0 MAXTAU] (default 50).
%
%    I = MUTUALINFORMATION(X, MAXTAU, FLAG) plots auto mutual information 
%    when FLAG is set to 1 (default 0).
%
%    I = MUTUALINFORMATION(X, MAXTAU, FLAG, BINS) uses the number BINS for
%    the estimation of the marginal probabilities. If not specified, the 
%    number of bins is computed using Freedman-Diaconis' rule.

% Copyright (c) 2020
% K. Hauke Kraemer, N. Marwan
% Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software and runs under MIT licence.

%% Assign input

y = varargin{1};
assert(isvector(y))

y = (y-mean(y)) ./ std(y);

try
    tau_max = varargin{2};
    assert(isscalar(tau_max), 'tau-values must be positive integers.')
    assert(tau_max > 0, 'tau-values must be positive integers.')
catch
    tau_max = 100;
end

try
    show = varargin{3};
    if ~(show==0 || show == 1)
        warning('input show needs to be 1 (display figure) or 0 (no figure displayed). Now set to 0.')
        show = 0;
    end
catch 
    show = 0;
end

try
    bins = varargin{4};
catch x 
    [~,edges] = histcounts(y);
    bins = length(edges)-1;
end

%% Check input
narginchk(1,4)
nargoutchk(0,1)

tau_min = 0;
Z(:,1)=tau_min:tau_max;

%% Compute mutual information
try
    I = mi(y, bins, tau_max,'nogui');
    Z(:,2) = zeros(1,size(I,3));
    Z(:,2) = I(1,1,:);
catch
    disp('no')
    Z(:,2) = NaN;
end

%% Plotting Automutualinformation against lag tau
if show
    figure
    plot(0:length(Z(:,2))-1, Z(:,2), '-.*', 'LineWidth', 2); hold on
    xlabel('time delay \tau')
    ylabel('mutual information [nats]')
    title('mutual information')
    set(gca, 'LineWidth', 2)
    grid on
end
