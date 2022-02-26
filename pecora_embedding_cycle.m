function [epsilon_mins, Y_old] = pecora_embedding_cycle(varargin)
% PECORA_EMBEDDING_CYCLE calculates the continuity statistics.
%    EPSILON_MINS = PECORA_EMBEDDING_CYCLE(X, TAU) computes vector EPSILON_MINS
%    containing the values of the continuity statistic based on Pecora et al., 
%    Chaos 17 (2007) for a given one- or multivariate time series provided by (N-by-M) 
%    matrix X, with N the length of the time series and M the number of time series, 
%    and for time delays given in vector TAU. In the first run, simply set TAU=0 to 
%    start with a non-embedded time series.
%
%    ... = PECORA_EMBEDDING_CYCLE(X, TAU, JS) computes the continuity statistic
%    using the time series from X identified by vector JS, which is a number from 1 to M
%    and according to the given delay in TAU.
%
%    ... = PECORA_EMBEDDING_CYCLE(X, TAU, JS, DELAY, SAMPLE_SIZE, THEILER, ALPHA, P, KNN, NORM)
%    computes the continuity statistic for possible delay times given by vector DELAY 
%    (default [0:50]), the fraction SAMPLE_SIZE of input time series length defining the 
%    size of the random phase space vector sample (default 0.1), the Theiler window THEILER 
%    (default 1), the significance level ALPHA (default 0.05) and binomial parameter P
%    (default 0.5) for the continuity statistic, the maximal number of points KNN in the 
%    delta-neighbourhood (default 13) which has to be at least 8, and the norm used for distance 
%    computations, which can be 'max' for Chebychev norm or 'euc' for Euclidean norm (default). 
%
%    [EPSILON_MINS, Y_OLD] = PECORA_EMBEDDING_CYCLE(...)
%    provides the phase space trajectory Y_OLD for which the additional tau's are being tested.
%
%    Further reading:
%    Pecora et al., Chaos 17 (2007)
%
%    See also UZAL_COST, PECUZAL_EMBEDDING

% Copyright (c) 2020
% K. Hauke Kraemer, N. Marwan
% Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software and runs under MIT licence.


%% in- and output check
narginchk(2,10)
nargoutchk(1,6)


%% Assign input

x = varargin{1};
% make the input time series a column vector
if size(x,1) < size(x,2)
    x = x';
end
% normalize time series
x = (x-mean(x)) ./ std(x);

taus = varargin{2};

try
    js = varargin{3};
catch
    js = ones(1, length(taus));
end

if length(taus) ~= length(js)
    error('tau inputs and js input must have equal length')
end

if max(js)>size(x,2)
    error('depicted time series in js do not match the time series input')
end
    

try
    delay_vals = varargin{4};
catch
    delay_vals = 0:50;
end
tN = length(delay_vals);


try
    sample_size = varargin{5};
catch
    sample_size = 0.1;
end
if sample_size < 0 || sample_size > 1
    warning('break percentage input must be a value in the interval [0 1]')
    sample_size = 0.1;
end

try
    theiler = varargin{6};
catch
    theiler = 1;
end

% confidence level for continuity statistic
try
    alpha = varargin{7};
catch
    alpha = 0.05;
end

% confidence level for continuity statistic
try
    p_val = varargin{8};
catch
    p_val = 0.5;
end

% delta neighborhood-sizes
try
    deltas = varargin{9};
catch
    deltas = 13;
end

methLib = {'euc', 'max'}; % the possible norms
try
    norm = varargin{10};
    if ~isa(norm,'char') || ~ismember(norm,methLib)
       warning(sprintf('Specified norm should be one of the following possible values: ''%s'', ''%s''.', methLib{:}))
       norm = 'euc';
    end
catch
    norm = 'euc';
end

%% Start computation

% Concerning the delta neighborhood of the continuity statistic:

% table for the continuity statistic
bino_table = get_binomial_table(p_val, alpha, deltas);
delta_points = bino_table(:,1);
epsilon_points = bino_table(:,2);

% considered neighbours
neighbours = delta_points(end); 

% intial phase space vector
for m = 1:length(taus)
    if m == 1
        Y_old = embed(x(:,js(m)), m, taus(m));
    else
        Y_old = embed2(Y_old, x(:,js(m)), taus(m));
    end
end
YN = length(Y_old);

% length of the reference point trajectory
NNN = floor(sample_size * (YN-delay_vals(end)));
% preallocate output
epsilon_mins = zeros(tN, size(x,2));

% select a random phase space vector sample
if sample_size == 1
    data_samps = 1:NNN;
else
    data_samps = datasample(1:YN-delay_vals(end), NNN, 'Replace', false); 
end

% loop over the different time series
for ts = 1:size(x, 2)
    % preallocate storing vector for continuity statistic
    epsilon_star = zeros(tN, NNN);
    
    % loop over all fiducial points
    for k = 1:NNN
        % bind the fiducial point from the trajectory sample
        fiducial_point = data_samps(k);
        
        % compute distances to all other points in dimension d.
        [distances, ~] = all_distances(Y_old(fiducial_point,:),...
                                      Y_old(1:end-delay_vals(end),:), norm);
        % sort these distances in ascending order
        [~,ind] = sort(distances);
        
        % temporal neighbours within Theiler window
        idx = max(1, fiducial_point - theiler):min(NNN, fiducial_point + theiler);
        % remove these neighbours from index list for distances
        ind(ismember(ind, idx)) = [];
        % construct neighbourhood
        NN_idxs = ind(1:neighbours)'; % indices of the valid neighbours 
        
        % preallocate storing vector for distances of epsilon neighborhood
        eps_distances = zeros(tN, neighbours);
        
        % loop over the different tau values
        for taus = 1:tN    
            tau = delay_vals(taus);
            eps_distances(taus,:) = abs(x(fiducial_point+tau,ts) - x(NN_idxs+tau,ts));
        end
       
        % now compute the minimum epsilon ranges for each delta 
        % neighbourhood size
        epsilon_star_delta = zeros(tN, length(delta_points));
        for i = 1:length(delta_points)
            l = delta_points(i);
            eps_range = sort(eps_distances(:,1:l),2);
            epsilon_star_delta(:,i) = eps_range(:, epsilon_points(i)); 
        end
        % Since we can not prefer one delta neighbourhood-size, we should 
        % take the minimum of all smallest scales.
        epsilon_star(:,k) = min(epsilon_star_delta, [], 2); 

    end
    % average over the fiducial points
    epsilon_min_avrg = mean(epsilon_star,2);

    % save all epsilon min vals corresponding to the different tau-vals for
    % this dimension-iteration
    epsilon_mins(:,ts) = epsilon_min_avrg;   
end
