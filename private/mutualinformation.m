function Z = mutualinformation(varargin)
%MUTUALINFORMATION computes the automutualinformation values of timeseries 
% `s` versus the time delayed version of itself.
%
% Minimum input-arguments : 1
% Maximum input-arguments : 4
%
%       MI = mutualinformation(s,tau_max,Show,bins)
% 
%
% Input arguments:
%
% s:        The time series, one wishes to compute the auto-
%           mutualinformation (MI) of.
% tau_max:  The maximum time lag up to which, starting from 0, MI gets
%           computed (Default is `tau_max` = 50).
% Show:     If set to 1, the results get displayed. (Default is `Show`= 0)
% bins:     Determines how many bins shall be used for the estimation of 
%           the marginal probabilities. If nothing is specified here the 
%           number of bins will be automatically computed using 
%           Freedman-Diaconis' rule.
%
% Output:
%
% MI:       The timelag-values (first column) and the corresponding 
%           mutualinformation-values (second column)
%
% This function uses mi()-function from N. Marwan's CRP toolbox.
%
% Copyright (c) 2020
% K. Hauke Kraemer, 
% Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software and runs under MIT licence.

%% Assign input

y = varargin{1};
assert(isvector(y))

y = (y-mean(y))/std(y);

try
    tau_max = varargin{2};
    assert(isscalar(tau_max),'tau-values must be positive integers.')
    assert(tau_max > 0,'tau-values must be positive integers.')
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
catch  
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
    I = mi(y',bins,tau_max);
    Z(:,2) = zeros(1,size(I,3));
    Z(:,2) = I(1,1,:);
catch
    disp('no')
    Z(:,2) = NaN;
end

%% Plotting Automutualinformation against lag tau
if show
    figure
    plot(0:length(Z(:,2))-1,Z(:,2),'-.*','LineWidth',2); hold on
    xlabel('time delay \tau')
    ylabel('mutual information [nats]')
    title('mutual information')
    set(gca,'LineWidth',2)
    set(gca,'FontSize',12)
    grid on
end

end

