function [epsilon_mins,gammas,dist_old,dist_,Y_old,fiducials] = pecora_embedding_cycle(varargin)
% PECORA_EMBED_TS is a unified approach to properly embed a time series based
% on the paper of Pecora et al., Chaos 17 (2007).
%
% Minimum input-arguments: 2
% Maximum input-arguments: 12
%
% [epsilon_mins, gammas, dist_old, dist, Y_old, fiducials] = 
%                       pecora_embedding_cycle(x,taus,js,delay_vals,...
%                                       datasample,theiler_window,...
%                                       undersampling,alpha,p_val,deltas,...
%                                       beta,norm);
%
% This function embeds the input time series `x` with different delay
% times tau. `x` can be a multivariate dataset, i.e. many univariate time 
% series in one matrix. The approach views the problem of choosing all 
% embedding parameters as being one and the same problem addressable using 
% a single statistical test formulated directly from the reconstruction 
% theorems. 
% This allows for varying time delays appropriate to the data and 
% simultaneously helps decide on embedding dimension. A second new 
% statistic, undersampling, acts as a check against overly long time delays 
% and overly large embedding dimension.
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
% `datasample`      Defines the size of the random phase space vector 
%                   sample the algorithm considers for each tau value, in 
%                   order to compute the continuity statistic. This is a
%                   float from the interval (0 1]. The size of the 
%                   considered sample is `datasample`*length of the current
%                   phase space trajectory (Default is 0.1, i.e. 1/10 of 
%                   the trajectory points will be considered).
% `theiler_window`  Defines a temporal correlation window for which no
%                   nearest neighbours are considered to be true, since
%                   they could be direct predecessors or successors of the
%                   fiducial point (Default is 1). When input `x` is a
%                   multivariate dataset, `theiler_window` needs to be
%                   chosen as the maximum theiler window from each of the
%                   time series.
% `undersampling`   A logical, set to 'true' the undersampling statistic
%                   also gets computed. Default is 'true'.
% `alpha`           The significance level for the continuity statistic.
% `p`               The binomial parameter used for obtaining the continuity 
%                   statistic (Default is `p=.5`).
% `deltas`          The maximal number of points in the delta-neighbourhood. 
%                   For a range of points from `8:deltas` the minimial 
%                   epsilon-scale according to the significance level `alpha` 
%                   and binomial parameter `p` gets computed and the minimum 
%                   of all considered delta-neighborhood-sizes is finally 
%                   taken. Must be at least 8 (Default is `deltas=8`)
% `beta`            The confidence level for the undersampling statistic.
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
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.

%% Assign input

% the input time series. Unlike the description in the docstring, yet there
% is just a univariate time series allowed
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

try 
    undersamp = varargin{7};
catch
    undersamp = true;
end

% confidence level for continuity statistic
try
    alpha = varargin{8};
catch
    alpha = 0.05;
end

% confidence level for continuity statistic
try
    p_val = varargin{9};
catch
    p_val = 0.5;
end

% delta neighborhood-sizes
try
    deltas = varargin{10};
catch
    deltas = 8;
end

% confidence level for undersampling statistic
try
    beta = varargin{11};
catch
    beta = 0.05;
end

methLib={'euc','max'}; % the possible norms
try
    norm = varargin{12};
    if ~isa(norm,'char') || ~ismember(norm,methLib)
       warning(['Specified norm should be one of the following possible values:',...
           10,sprintf('''%s'' ',methLib{:})])
       norm = 'euc';
    end
catch
    norm = 'euc';
end


% Matlab in- and output check
narginchk(2,12)
nargoutchk(1,6)

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


% select a random phase space vector sample. One could of course 
% take all points from the actual trajectory here, but this would
% be computationally overwhelming. This is why I incorporated the
% input parameter 'sample_size' as the fraction of all phase space
% points.
if sample_size == 1
    data_samps = 1:NNN;
else
    data_samps = datasample(1:YN-delay_vals(end),...
        NNN,'Replace',false);
    
end
fiducials = data_samps;

% loop over the different time series
for ts = 1:size(x,2)
    
    % preallocate storing vector for minimum epsilon from the
    % continuity statistic
    epsilon_star = zeros(tN,NNN);
    
    if undersamp
        % for the undersampling statistic perform convolution of the time
        % series
        [xx,sigma] = convolute(x(:,ts),x(:,ts));
        
        % preallocate storing vector for maximum gamma from the
        % undersampling statisctic
        gamma_k = zeros(tN,NNN);
        
        dist = NaN*ones(tN,NNN);
    end
    
    % loop over all fiducial points, look at their neighbourhoods and
    % compute the continuity as well as the undersampling statistic for
    % each of them. Eventually we average over all the values to get
    % the continuity and undersampling statistic for the actual tau
    % value.
    for k = 1:NNN

        % bind the fiducial point from the trajectory sample
        fiducial_point = data_samps(k);
        

        % compute distances to all other points in dimension d. You'll
        % find the helper function at the end of the script
        [distances, ~] = all_distances(Y_old(fiducial_point,:),...
                                                Y_old(1:end-delay_vals(end),:),norm);
                                       
        % sort these distances in ascending order
        [~,ind] = sort(distances);
        
        % 1) compute the continuity statistic

        % loop over all delta-neighbourhoods. Here we take the table
        % from the paper, which is based on an confidence level alpha =
        % 0.05. In the final implementation one should be able to
        % choose an alpha, maybe at least from 0.05 or 0.01. Therefore
        % one could look up a similar table (binomial distribution). We
        % try all these number of points, i.e. we try many different
        % delta's, as mentioned in the paper. For each of the delta
        % neighbourhoods we loop over the epsilons (decreasing values)
        % and stop, when we can not reject the null anymore. After
        % trying all deltas we pick the maximum epsilon from all
        % delta-trials. This is then the final epsilon for one specific
        % tau and one specific fiducial point. Afterwards we average
        % over all fiducial points.  
        
            
        % set index-counter for the upcoming for-loop over the different tau 
        % values
        tau_counter = 1;
        
        % preallocate storing vector for distances of epsilon neighborhood
        eps_distances = zeros(tN,neighbours);
        
        % loop over the different tau values. Starting at 0 is important,
        % especially when considering a mutlivariate input data set to choose
        % from (in the final implemenation)
        
        for taus = 1:tN
            
            tau = delay_vals(taus);
            
            % create new phase space vector. 'embed2' is a helper function, 
            % which you find at the end of the script. At this point one could
            % incorporate onther loop over all the input time series, when
            % allowing multivariate input. 
            Y_new = embed2(Y_old,x(:,ts),tau);
                      
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

                    % save the last component of this 
                    % neighbour for later comparison to the epsilon
                    % interval. Therefore, first check that the 
                    % neighbour is not in the temporal correlation 
                    % window (determined by the input 'theiler').
                    % If it is, go on to the next.
                    if ind(l) > fiducial_point + theiler || ...
                            ind(l) < fiducial_point - theiler
                        % save the distance of the valid neighbour 
                        eps_distances(tau_counter,nei) = ...
                            abs(Y_new(ind(l),size(Y_new,2))...
                            -Y_new(fiducial_point,size(Y_new,2)));
                        if nei == 1
                            dist_old(ts,k)=distances(ind(l));
                        end
                        break
                    else
                        % check the next neighbour
                        l = l + 1;
                    end
                    % make sure the data set is sufficiently
                    % sampled, if not pass an error. Since 'ind' is
                    % a vector storing all indices of nearest
                    % neighbours and l is exceeding its length
                    % would mean that all neighbours are so close
                    % in time that they fall within the correlation
                    % window OR the number of neighbours one is
                    % looking at (see table 1 in the paper) is too
                    % large
                    if l > length(ind)
                        error('not enough neighbours')
                    end
                end 
                
                % 2) compute undersampling statistic
                if nei == 1 && undersamp
                    % undersampling statistic

                    % compute new neighbourhood
                    [distances2,comp_dist_new] = ...
                        all_distances(Y_new(fiducial_point,:),Y_new,norm);

                    [~,ind2] = sort(distances2);

                    % search for valid neighbour w.r.t. theiler window
                    for i = 1:length(ind2)
                        if ind2(i) > fiducial_point + theiler || ...
                            ind2(i) < fiducial_point - theiler
                            
                            % bind the distance
                            dist(tau_counter,k) = distances2(ind2(i));                            
                            break
                        else
                            continue
                        end
                    end

                    if isnan(dist(tau_counter,k))
                        error('no valid neighbours, please reduce the theiler wondow or increase data size')
                    end

                    % perform undersampling statistic test. You'll find the 
                    % helper function undersampling at the end of this script.
                  
                    % now run the undersampling statistic function on the
                    % already computed convolution sigma and for the
                    % distances stored in dist under the significance level
                    % alpha
                    if length(taus) == 1
                        gg = undersampling(xx,sigma,comp_dist_new(ind2(i),end-1:end),beta);
                        gamma_k(tau_counter,k) = max(gg);
                    else
                        gamma_k(tau_counter,k) = undersampling(xx,sigma,comp_dist_new(ind2(i),end),beta);
                    end
                end
                
                % check the next neighbour
                l = l + 1;
                   
            end 
                
            % increase counter for tau-loop
            tau_counter = tau_counter + 1;  
            
        end
        
        % now compute the minimum epsilon ranges for each delta 
        % neighbourhood size.

        % if the number of neighbours from the delta
        % neighbourhood, which get projected into the epsilon
        % neighbourhood (Fig.1 in the paper) are smaller than
        % the amount needed to reject the null, we break and
        % store this particular epsilon value (which 
        % determines the epsilon neighbourhood-size)

        epsilon_star_delta = zeros(tN,length(delta_points));
        for i = 1:length(delta_points)
            l = delta_points(i);
            eps_range = sort(eps_distances(:,1:l),2);   
            epsilon_star_delta(:,i) = eps_range(:,epsilon_points(i));         
            clear eps_range
        end


        % In the paper it is not clearly stated how to preceed here. We
        % are looking for the smallest scale for which we can not
        % reject the null. Since we can not prefer one delta
        % neighbourhood-size, we should take the maximum of all
        % smallest scales.
        epsilon_star(:,k) = min(epsilon_star_delta,[],2); 
                                  
    end
    
    if ~undersamp
        dist_{ts} = NaN;
    else
        dist_{ts} = dist;
    end
    
    clear dist
    
    % average over the fiducial points
    % 1) continuity statistic
    epsilon_min_avrg = mean(epsilon_star,2);
    
    % save all epsilon min vals corresponding to the different tau-vals for
    % this dimension-iteration
    epsilon_mins(:,ts) = epsilon_min_avrg;
    
    % 2) gamma statistic
    if undersamp
        
        gamma_avrg = mean(gamma_k,2);
        % save all gamma vals corresponding to the different tau-vals for
        % this dimension-iteration
        gammas(:,ts) = gamma_avrg;
    end
       
end
        
end