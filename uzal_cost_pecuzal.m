function L_decrease = uzal_cost_pecuzal(Y, Y_trial, varargin)
% UZAL_COST_PECUZAL Calculates the maximum L-decrease of two trajectories.
%    L_decrease = UZAL_COST_PECUZAL(Y1, Y2) computes vector L_decrease 
%    containing the maximum decrease of the cost function L proposed in 
%    Uzal et al., Phys. Rev. E 84, 016223, 2011 of a given phase space 
%    trajectory provided by (N-by-M) matrix Y1 and a given phase space 
%    trajectory provided by (N'-by-M') matrix Y2 with N, N' the length
%    of phase space trajectories and M, M' its dimensions.
%
%    L_decrease = UZAL_COST_PECUZAL(Y1, Y2, TW) computes vector L_decrease 
%    containing the maximum decrease of the cost function L of Y1 and Y2 
%    considering maximal time horizons up to TW (default 50)
%
%    L_decrease = UZAL_COST_PECUZAL(..., Name, Value) specifies further optional 
%    parameters for the algorithm using one or more Name, Value pair arguments.
%
%    Optional name-value-arguments:
%      'theiler'     - (default 1) the temporal correlation window (Theiler window) 
%                      for which no nearest neighbours are considered, because they 
%                      could be direct predecessors or successors of the fiducial point. 
%                      When input X is a multivariate dataset, the Theiler window needs 
%                      to be chosen as the maximum from each of the time series.
%      'k'           - (default 3) considered number of nearest neighbors for Uzal's
%                      L-statistic.
%      'econ'        - (default False) Economy-mode for L-statistic computation. Instead of
%                      computing L-statistics for time horizons `2:Tw`, here we only compute them for
%                      `2:2:Tw`.
%
%    Further reading:
%    Uzal et al., Phys. Rev. E 84, 016223, 2011
%
%    See also PECUZAL_EMBEDDING, PECORA_EMBEDDING_CYCLE

% Copyright (c) 2020
% K. Hauke Kraemer, N. Marwan
% Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software and runs under MIT licence.

%% in- and output check
narginchk(1,9)
nargoutchk(1,1)

%% Assign input

% default values
Tw = 50;
theiler = 1;
k = 3;
econ = false;

% required and optional arguments
p = inputParser;

validScalarPosNum1 = @(x) isnumeric(x) && isscalar(x) && (x > 0) && rem(x,1)==0;
validScalarPosNum2 = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validDimension = @(x) isnumeric(x) && ismatrix(x);
validType = @(x) islogical(x);

addRequired(p,'Y',validDimension);
addRequired(p,'Y_trial',validDimension);
addOptional(p,'Tw',Tw,validScalarPosNum2);
addParameter(p,'theiler',theiler,validScalarPosNum1);
addParameter(p,'k',k,validScalarPosNum2);
addParameter(p,'econ',econ,validType);

% parse input arguments
parse(p, Y, Y_trial, varargin{:})

% assign variables with the resulting argument input
Y = p.Results.Y;
Y_trial = p.Results.Y_trial;
Tw = p.Results.Tw;
theiler = p.Results.theiler;
k = p.Results.k;
econ = p.Results.econ;

norm = 'euc';
metric = 'euclidean';

% make the input time series row vectors
if size(Y, 1) < size(Y, 2)
    Y = Y';
end

if size(Y_trial, 1) < size(Y_trial, 2)
    Y_trial = Y_trial';
end

%% Start computation

% select a random phase space vector sample. 
NNN = length(Y_trial)-1;
D1 = size(Y,2);
D2 = size(Y_trial,2);

if econ
    tws = 2:2:Tw; % start at 2 will eliminate bad results for noise
else
    tws = 2:Tw; % start at 2 will eliminate bad results for noise
end

% preallocation
epsilon_k2 = zeros(1,NNN);          % neighborhood size
epsilon_k2_trial = zeros(1,NNN);    % neighborhood size
E2 = zeros(NNN,length(tws));                 % conditional variance
E2_trial = zeros(NNN,length(tws));           % conditional variance
neighborhood = zeros(k+1,D1);       % epsilon neighbourhood
neighborhood_trial = zeros(k+1,D2); % epsilon neighbourhood

% initial L-decrease
dist_former = 99999999;

% loop over each time horizon
cnt = 1;
for j = 1:length(tws) 
    T = tws(j);
    NN = length(Y_trial)-T;
    if NN < 1
        error("Time series too short for given possible delays and Theiler window to find enough nearest neighbours")
    end
    
    for fiducial_point = 1:NN

        % find nearest neighbours for the fiducial point and
        % compute distances to all other points. 
        distances = all_distances(Y(fiducial_point,:), Y(1:NN,:), norm);                                      
        % sort these distances in ascending order
        [~,ind] = sort(distances);
        % temporal neighbours within Theiler window
        idx = max(1,fiducial_point - theiler):min(NN,fiducial_point + theiler);
        % remove these neighbours from index list for distances
        ind(ismember(ind, idx)) = [];
        
        distances_trial = all_distances(Y_trial(fiducial_point,:), Y_trial(1:NN,:), norm);                                      
        % sort these distances in ascending order
        [~,ind_trial] = sort(distances_trial);
        % remove these neighbours from index list for distances
        ind_trial(ismember(ind_trial, idx)) = [];
        
        %% construct neighbourhoods
        neighborhood(1,:) = Y(fiducial_point,:);
        neighborhood(2:k+1,:) = Y(ind(1:k),:); % values of the neighbour vectors
        neighborhood_trial(1,:) = Y_trial(fiducial_point,:);
        neighborhood_trial(2:k+1,:) = Y_trial(ind_trial(1:k),:); % values of the neighbour vectors
        eps_idx(1) = fiducial_point;
        eps_idx(2:k+1) = ind(1:k)'; % indices of the valid neighbours 
        eps_idx_trial(1) = fiducial_point;
        eps_idx_trial(2:k+1) = ind_trial(1:k)'; % indices of the valid neighbours
   
        %% estimate size of neighbourhood
        pd = pdist(neighborhood, metric);
        epsilon_k2(fiducial_point) = (2/(k*(k+1))) * sum(pd.^2);  % Eq. 16
        pd_trial = pdist(neighborhood_trial, metric);
        epsilon_k2_trial(fiducial_point) = (2/(k*(k+1))) * sum(pd_trial.^2);  % Eq. 16
    
        %% estimate E2
        % determine neighborhood T timesteps ahead
        eps_ball = Y(eps_idx+T,:);
        eps_ball_trial = Y_trial(eps_idx_trial+T,:);
        % compute center of mass
        u_k = sum(eps_ball) / (k+1);  % Eq.14
        u_k_trial = sum(eps_ball_trial) / (k+1);  % Eq.14
        % compute E_k2
        E2(fiducial_point,cnt) = sum(rms(eps_ball-u_k).^2); % Eq.13 
        E2_trial(fiducial_point,cnt) = sum(rms(eps_ball_trial-u_k_trial).^2); % Eq.13 

    end
    % compute distance of L-values and check whether that distance can be
    % increased
    
    % Average E2 over all fiducial points         
    E2_avrg = mean(E2(1:NN,1:cnt),2);                   % Eq. 15
    E2_avrg_trial = mean(E2_trial(1:NN,1:cnt),2);
    sigma2 = E2_avrg ./ epsilon_k2(1:NN)'; % noise amplification σ², Eq. 17
    sigma2_trial = E2_avrg_trial ./ epsilon_k2_trial(1:NN)'; 
    sigma2_avrg = mean(sigma2); % averaged value of the noise amplification, Eq. 18
    sigma2_avrg_trial = mean(sigma2_trial);
    alpha2 = 1 / mean(epsilon_k2(1:NN).^(-1)); % for normalization, Eq. 21
    alpha2_trial = 1 / mean(epsilon_k2_trial(1:NN).^(-1));
    L = log10(sqrt(sigma2_avrg)*sqrt(alpha2));
    L_trial = log10(sqrt(sigma2_avrg_trial)*sqrt(alpha2_trial));
    
    dist = L_trial - L;
    if isnan(dist)
        error("Computed 0-distances. You might use model-data, thus try to add minimal additive noise to the signal you wish to embed and try again.")
    end
    if (dist > dist_former) && (dist_former<0)
        break
    else
        dist_former = dist;
    end
    cnt = cnt + 1;
end

L_decrease = dist_former;

end
