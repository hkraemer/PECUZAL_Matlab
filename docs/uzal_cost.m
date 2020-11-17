function L = uzal_cost(varargin)
% UZAL_COST computes the cost function proposed in Uzal et al., Phys. Rev.
% E 84, 016223, 2011 of a given Phase space trajectory.
%
% Minimum input: 1
% Maximum input: 6
% L = uzal_cost(Y,theiler,k,Tw,sample_size,norm)
%
%
% Input:
%
% Y:                  state space vector
% theiler:            Theiler window for excluding serial correlated points 
%                     from neighbourhood (Default is theiler=1)
% k:                  number of nearest neighbors to be considered (Default
%                     is k=3
% Tw:                 Time forward parameter (Default is Tw = 40)
% sample_size:        Number of considered fiducial points as a fraction of 
%                     input time series length (Default is sample_size = 0.2)
% norm:               The norm used for distance computations. Choose from
%                     'max' (Chebychev norm) and 'euc' (Euclidean norm).
%                     Default is norm = 'euc'
%
%
% Output:
%
% L:                  The value of the proposed cost-function
%
%
% Copyright (c) 2020
% K. Hauke Kraemer, 
% Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software and runs under MIT licence.

%% Assign input

% the input time series. Unlike the description in the docstring, yet there
% is just a univariate time series allowed
Y = varargin{1};

% make the input time series row vectors
if size(Y,1)<size(Y,2)
    Y = Y';
end

try
    theiler = varargin{2};
catch
    theiler = 1;
end

try
    k = varargin{3};
catch
    k = 3;
end

try
    Tw = varargin{4};
catch
    Tw = 40;
end


try
    sample_size = varargin{5};
catch
    sample_size = 0.2;
end
if sample_size < 0 || sample_size > 1
    warning('sample size input must be a value in the interval [0 1]')
    sample_size = 0.2;
end

methLib={'euc','max'}; % the possible norms
try
    norm = varargin{6};
    if ~isa(norm,'char') || ~ismember(norm,methLib)
       warning(['Specified norm should be one of the following possible values:',...
           10,sprintf('''%s'' ',methLib{:})])
    end
catch
    norm = 'euc';
end

if strcmp(norm,'euc')
    metric = 'euclidean';
elseif strcmp(norm,'max')
    metric = 'chebychev';
end

% Matlab in- and output check
narginchk(1,6)
nargoutchk(1,1)

%% Start computation

% select a random phase space vector sample. 
NN = length(Y)-Tw;
NNN = floor(sample_size*NN);
data_samps = datasample(1:NN,NNN,'Replace',false);

% preallocation
E_k2_avrg = zeros(1,NNN);
epsilon_k2 = zeros(1,NNN);
% loop over all fiducial points
for ks = 1:NNN
    % bind the fiducial point from the trajectory sample
    fiducial_point = data_samps(ks);
    
    % find nearest neighbours for the fiducial point and
    % compute distances to all other points in dimension d. 
    [distances, ~] = all_distances(Y(fiducial_point,:),Y(1:end-Tw,:),norm);
                                       
    % sort these distances in ascending order
    [~,ind] = sort(distances);
    
    %% construct neighbourhood
    % loop over all neighbours and get their distances to the
    % fiducial point and neighbor-indices
    eps_idx = zeros(1,k+1);
    eps_idx(1) = fiducial_point;
    neighborhood = zeros(k+1,size(Y,2));
    neighborhood(1,:) = Y(fiducial_point,:);
 
    l = 2;  % start with the first neighbour which is not 
            % the fiducial point itself
    for nei = 1:k
        % this while loop gurantees, that we look at a true
        % neighbour and not a one which lies in the
        % correlation window of the fiducial point
        flag = true;
        while flag
            if ind(l) > fiducial_point + theiler || ...
                    ind(l) < fiducial_point - theiler
                % save the index of the valid neighbour 
                eps_idx(nei+1) = ind(l);
                % save vector to neighborhood
                neighborhood(nei+1,:) = Y(ind(l),:);
                l = l + 1;
                flag = false;
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
    
    %% estimate size of neighbourhood
    pd = pdist(neighborhood,metric);
    epsilon_k2(ks) = (2/(k*(k+1))) * sum(pd.^2);  % Eq. 16
    
    %% estimate E_k2
    % loop over the different Time horizons
    E_k2 = zeros(1,Tw);
    for T = 1:Tw
        % determine neighborhood T timesteps ahead
        eps_ball = Y(eps_idx+T,:);
        % compute center of mass
        u_k = sum(eps_ball) / (k+1);  % Eq.14
        % compute E_k2
        if strcmp(norm,'euc')
            E_k2(T) = sum((sqrt(sum((eps_ball-u_k).^2))).^2) / (k+1);  % Eq.13 
        elseif strcmp(norm,'max')
            E_k2(T) = sum((max(abs(eps_ball-u_k))).^2) / (k+1);  % Eq.13 
        end 
    end
    % time average over all prediction horizons
    E_k2_avrg(ks) = mean(E_k2); % Eq. 15       
end
% compute sigma_k2
sigma_k2 = E_k2_avrg ./ epsilon_k2; % Eq. 17

sigma_k2_avrg = mean(sigma_k2);   % Eq. 18

alpha_k2 = 1/sum(1./epsilon_k2); % Eq. 21

% Output
L = log10(sqrt(sigma_k2_avrg)*sqrt(alpha_k2));  % Eq. 26

end
