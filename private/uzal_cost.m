function L = uzal_cost(varargin)
% UZAL_COST Calculates the values of L-statistics.
%    L = UZAL_COST(X) computes vector L containing the values of the cost 
%    function proposed in Uzal et al., Phys. Rev. E 84, 016223, 2011 of a 
%    given phase space trajectory provided by (N-by-M) matrix X, with N the length
%    of phase space trajectory and M the dimension.
%
%    L = UZAL_COST(X, THEILER, K, TW, SAMPLE_SIZE, NORM) calculates L using
%    Theiler window THEILER (default 1), for excluding serial correlated points,
%    the number K of nearest neighbors to be considered (default 3), the time
%    forward parameter TW (default 40), the fraction SAMPLE_SIZE of input time 
%    series length defining the number of considered fiducial points (default 0.2),
%    and the norm used for distance computations, which can be 'max' for Chebychev norm
%    or 'euc' for Euclidean norm (default). 

% Copyright (c) 2020
% K. Hauke Kraemer, 
% Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software and runs under MIT licence.

%% in- and output check
narginchk(1, 6)
nargoutchk(1, 1)


%% Assign input

% the input time series. Unlike the description in the docstring, yet there
% is just a univariate time series allowed
Y = varargin{1};

% make the input time series row vectors
if size(Y, 1) < size(Y, 2)
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
    if ~isa(norm, 'char') || ~ismember(norm, methLib)
       warning(sprintf('Specified norm should be one of the following possible values: ''%s'', ''%s''.', methLib{:}))
    end
catch
    norm = 'euc';
end

if strcmp(norm,'euc')
    metric = 'euclidean';
elseif strcmp(norm,'max')
    metric = 'chebychev';
end

%% Start computation

% select a random phase space vector sample. 
NN = length(Y)-Tw;
NNN = floor(sample_size*NN);
data_samps = datasample(1:NN, NNN, 'Replace', false);

use_vectorized = 0;
if NNN > 2500 
   use_vectorized = 1;
end

% preallocation
E_k2_avrg = zeros(1,NNN);
epsilon_k2 = zeros(1,NNN);
% loop over all fiducial points
for ks = 1:NNN
    % bind the fiducial point from the trajectory sample
    fiducial_point = data_samps(ks);
    
    % find nearest neighbours for the fiducial point and
    % compute distances to all other points in dimension d. 
    distances = all_distances(Y(fiducial_point,:), Y(1:end-Tw,:), norm);
                                       
    % sort these distances in ascending order
    [~,ind] = sort(distances);
 
    % temporal neighbours within Theiler window
    idx = max(1,fiducial_point - theiler):min(NN,fiducial_point + theiler);
    % remove these neighbours from index list for distances
    ind(ismember(ind, idx)) = [];

    %% construct neighbourhood
    neighborhood(1,:) = Y(fiducial_point,:);
    neighborhood(2:k+1,:) = Y(ind(1:k),:); % values of the neighbour vectors
    eps_idx(1) = fiducial_point;
    eps_idx(2:k+1) = ind(1:k)'; % indices of the valid neighbours 
   
    %% estimate size of neighbourhood
    pd = pdist(neighborhood,metric);
    epsilon_k2(ks) = (2/(k*(k+1))) * sum(pd.^2);  % Eq. 16
    
    %% estimate E_k2
    % loop over the different Time horizons
    E_k2 = zeros(1,Tw);
    
    % series of neighbours
    idx = [eps_idx' ones(k+1,1)]*[ones(1,Tw); 1:Tw]; 

    if use_vectorized
        % neighborhood T timesteps ahead
        eps_ball = reshape(Y(idx,:),k+1, Tw,size(Y,2));
        eps_ball = permute(eps_ball,[1 3 2]);
        % center of mass
        u_k=(sum(eps_ball,1)) / (k+1);
        % compute E_k2
        if strcmp(norm,'euc')
            E_k2 = squeeze(sum(sqrt(sum((eps_ball - u_k).^2)).^2) / (k+1));
        elseif strcmp(norm,'max')
            E_k2 = squeeze(sum(max(abs(eps_ball - u_k)).^2) / (k+1));
        end

    else
        for T = 1:Tw
            % determine neighborhood T timesteps ahead
            eps_ball = Y(eps_idx+T,:);
            % compute center of mass
            u_k = sum(eps_ball) / (k+1);  % Eq.14
            % compute E_k2
            if strcmp(norm,'euc')
                E_k2(T) = sum(sqrt(sum((eps_ball-u_k).^2)).^2) / (k+1);  % Eq.13 
            elseif strcmp(norm,'max')
                E_k2(T) = sum(max(abs(eps_ball-u_k)).^2) / (k+1);  % Eq.13 
            end 
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
