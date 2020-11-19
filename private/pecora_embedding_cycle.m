function [epsilon_mins,gammas,dist_old,dist_,Y_old,fiducials] = pecora_embedding_cycle(varargin)
% PECORA_EMBEDDING_CYCLE calculates the continuity statistics.
%    PECORA_EMBEDDING_CYCLE computes the continuity statistic `epsilon_mins`
% based on the paper of Pecora et al., Chaos 17 (2007).
%
% Minimum input-arguments: 2
% Maximum input-arguments: 12
%
% [epsilon_mins, gammas, dist_old, dist, Y_old, fiducials] = 
%                       pecora_embedding_cycle(x,taus,js,delay_vals,...
%                                       sample_size,theiler,alpha,p,...
%                                       KNN,norm);
%
% Input:    
%
% `x`               A uni- or multivariate time series, which needs to be
%                   embedded. If the input data is a multivariate set, the
%                   algorithm scans all time series and constructs the
%                   according statistics.
% `taus`            An array, which stores the taus for embedding the time
%                   series in `x`. In the first run, simply set `taus=[0]`,
%                   to start with a non-embedded time series.
% `js`              An array, which stores the time series numbers as given
%                   in 'x' according to the given 'taus'.
% `delay_vals`      A vector defining for which possible delay times tau the
%                   algorithm shall look (Default is `delay_vals` = 0:50).
% `sample_size`     Defines the size of the random phase space vector 
%                   sample the algorithm considers for each tau value, in 
%                   order to compute the continuity statistic. This is a
%                   float from the interval (0 1]. The size of the 
%                   considered sample is `sample_size`*length of the current
%                   phase space trajectory (Default is 0.1, i.e. 1/10 of 
%                   the trajectory points will be considered).
% `theiler`         Defines a temporal correlation window for which no
%                   nearest neighbours are considered to be true, since
%                   they could be direct predecessors or successors of the
%                   fiducial point (Default is 1). When input `x` is a
%                   multivariate dataset, `theiler` needs to be chosen as 
%                   the maximum theiler window from each of the time series..
% `alpha`           The significance level for the continuity statistic.
% `p`               The binomial parameter used for obtaining the continuity 
%                   statistic (Default is `p=.5`).
% `KNN`             The maximal number of points in the delta-neighbourhood. 
%                   For a range of points from `8:KNN` the minimial 
%                   epsilon-scale according to the significance level `alpha` 
%                   and binomial parameter `p` gets computed and the minimum 
%                   of all considered delta-neighborhood-sizes is finally 
%                   taken. Must be at least 8 (Default is `KNN=8`)
% `norm`            norm for distance calculations. 'euc' (Default) or 'max' 
%
% Output: 
%
% `epsilon_mins`    Continuity statistic. A matrix storing all epsilon-
%                   mins as a function of `delay_vals` for each encountered 
%                   time series in input `x`.
% `gammas`          Undersampling statistic. If `undersampling` input is 
%                   'true' this statistic will be returned. A matrix 
%                   storing all gamma values as a function of `delay_vals` 
%                   for each encountered time series in input `x`.
% `dist_old`        All nearest neighbour distances to each fiducial point
%                   for the given trajectory. This is a number of fiducial 
%                   points - dimensional array, for each encountered time
%                   series.
% `dist`            All nearest neighbour distances to each fiducial point
%                   for each value of considered tau. This is a 
%                   `length(delay_vals)`*number of fiducial points - 
%                   dimensional array for each encountered time series. 
%                   Each of them is stored as entry in a cell array.
% `Y_old`           The phase space trajectory for which the additional
%                   tau's are being tested
% `fiducials`       The indices of the chosen fiducial points


% Copyright (c) 2020
% K. Hauke Kraemer, 
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
if size(x,1)<size(x,2)
    x = x';
end
% normalize time series
x = (x-mean(x))./std(x);

taus = varargin{2};

try
    js = varargin{3};
catch
    js = ones(1,length(taus));
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
    deltas = 8;
end

methLib={'euc','max'}; % the possible norms
try
    norm = varargin{10};
    if ~isa(norm,'char') || ~ismember(norm,methLib)
       warning(['Specified norm should be one of the following possible values:',...
           10,sprintf('''%s'' ',methLib{:})])
       norm = 'euc';
    end
catch
    norm = 'euc';
end

%% Start computation

% Concerning the delta neighborhood of the continuity statistic:

% table for the continuity statistic
bino_table = get_binomial_table(p_val,alpha,deltas);
delta_points = bino_table(:,1);
epsilon_points = bino_table(:,2);

% considered neighbours
neighbours = delta_points(end); 

% intial phase space vector
for m = 1:length(taus)
    if m == 1
        Y_old = embed(x(:,js(m)),m,taus(m));
    else
        Y_old = embed2(Y_old,x(:,js(m)),taus(m));
    end
end
YN = length(Y_old);
% length of the reference point trajectory
NNN = floor(sample_size*(YN-delay_vals(end)));

% preallocate output
epsilon_mins = zeros(tN,size(x,2));
gammas = zeros(tN,size(x,2));
dist_ = cell(1,size(x,2));
dist_old = NaN*ones(size(x,2),NNN);

% select a random phase space vector sample.
if sample_size == 1
    data_samps = 1:NNN;
else
    data_samps = datasample(1:YN-delay_vals(end),...
        NNN,'Replace',false); 
end
fiducials = data_samps;

% loop over the different time series
for ts = 1:size(x,2)
    
    % preallocate storing vector for continuity statistic
    epsilon_star = zeros(tN,NNN);

    % preallocate nearest neighbour indices matrix
    NN_idxs = zeros(NNN,neighbours);
    
    % loop over all fiducial points
    for k = 1:NNN
        % bind the fiducial point from the trajectory sample
        fiducial_point = data_samps(k);
        
        % compute distances to all other points in dimension d.
        [distances, ~] = all_distances(Y_old(fiducial_point,:),...
                                      Y_old(1:end-delay_vals(end),:),norm);                                      
        % sort these distances in ascending order
        [~,ind] = sort(distances);
        
        % set index-counter for the upcoming for-loop over the different tau 
        % values
        tau_counter = 1;
        
        % preallocate storing vector for distances of epsilon neighborhood
        eps_distances = zeros(tN,neighbours);
        
        % loop over the different tau values
        for taus = 1:tN    
            tau = delay_vals(taus);
            % create new phase space vector.
            Y_new = embed2(Y_old,x(:,ts),tau);          
            if taus == 1
                % loop over all neighbours (determined by the delta-
                % neighbourhood-size) and get their distances to the
                % fiducial point. 
                l = 2;  % start with the first neighbour which is not 
                        % the fiducial point itself
                for nei = 1:neighbours
                    % this while loop gurantees, that we look at a true
                    % neighbour and not a one which lies in the
                    % correlation window of the fiducial point
                    while true
                        if ind(l) > fiducial_point + theiler || ...
                                ind(l) < fiducial_point - theiler
                             NN_idxs(k,nei) = ind(l);
                             l = l + 1;
                             break
                        else
                            % check the next neighbour
                            l = l + 1;
                        end
                        % make sure the data set is sufficiently
                        % sampled
                        if l > length(ind)
                            error('not enough neighbours')
                        end
                    end 
                end
            end 
            
            % save the distance of the valid neighbour 
            eps_distances(tau_counter,:) = ...
                abs(Y_new(NN_idxs(k,:),size(Y_new,2))...
                -Y_new(fiducial_point,size(Y_new,2)));
            % increase counter for tau-loop
            tau_counter = tau_counter + 1;  
        end
       
        % now compute the minimum epsilon ranges for each delta 
        % neighbourhood size.
        epsilon_star_delta = zeros(tN,length(delta_points));
        for i = 1:length(delta_points)
            l = delta_points(i);
            eps_range = sort(eps_distances(:,1:l),2);   
            epsilon_star_delta(:,i) = eps_range(:,epsilon_points(i));         
            clear eps_range
        end
        % Since we can not prefer one delta neighbourhood-size, we should 
        % take the minimum of all smallest scales.
        epsilon_star(:,k) = min(epsilon_star_delta,[],2);                                
    end
    clear dist
    
    % average over the fiducial points
    epsilon_min_avrg = mean(epsilon_star,2);
    % save all epsilon min vals corresponding to the different tau-vals for
    % this dimension-iteration
    epsilon_mins(:,ts) = epsilon_min_avrg;   
end
        
end
