function [Y_tot,tau_vals,ts_vals,epsilons,LS] = pecuzal_embedding(x,varargin)
% PECUZAL_EMBEDDING is fully automated phase space reconstruction method
% based on the ideas of Pecora et al., Chaos 17 (2007) & Uzal et al., Phys.
% Rev. E (2011).
%
% Minimum input-arguments: 1
% Maximum input-arguments: 21
%
% [Y, tau_vals, ts_vals, epsilon_mins, Ls] = 
%                       pecuzal_embedding(x, tau_vals [,Name, Value]);
%
% This function embeds the input time series `x` with different delay
% times tau. `x` can be a multivariate dataset, i.e. many univariate time 
% series in one matrix. The approach views the problem of choosing all 
% embedding parameters as being one and the same problem addressable using 
% a single statistical test formulated directly from the reconstruction 
% theorems. 
% In combination with the L-statistic, a cost-function, which evaluates the
% goodness of a reconstruction, this allows for varying time delays 
% appropriate to the data and simultaneously helps decide on embedding 
% dimension. 
%
% Input:    
% 
% Required:
% `x`               A uni- or multivariate time series, which needs to be
%                   embedded. If the input data is a multivariate set, the
%                   algorithm scans all time series and constructs the
%                   according statistics.
% Optional:
% `delay_vals`      A vector defining for which possible delay times tau the
%                   algorithm shall look (Default is `delay_vals` = 0:50).
%
% Optional name-value-arguments:
% `datasample`      Defines the size of the random phase space vector 
%                   sample the algorithm considers for each tau value, in 
%                   order to compute the continuity statistic. This is a
%                   float from the intervall (0 1]. The size of the 
%                   considered sample is `datasample`*length of the current
%                   phase space trajectory (Default is 1, i.e. all 
%                   the trajectory points will be considered).
% `theiler`         Defines a temporal correlation window for which no
%                   nearest neighbours are considered to be true, since
%                   they could be direct predecessors or successors of the
%                   fiducial point (Default is 1). When input `x` is a
%                   multivariate dataset, `theiler` needs to be
%                   chosen as the maximum `theiler` from each of the
%                   time series.
% `alpha`           The significance level for the continuity statistic
%                   (Default is `alpha=0.05`)
% `p`               The binomial parameter used for obtaining the continuity 
%                   statistic (Default is `p=.5`).
% `deltas`          The maximal number of points in the delta-neighbourhood. 
%                   For a range of points from `8:deltas` the minimial 
%                   epsilon-scale according to the significance level `alpha` 
%                   and binomial parameter `p` gets computed and the minimum 
%                   of all considered delta-neighborhood-sizes is finally 
%                   taken. Must be at least 8 (Default is `deltas=13`) 
% `k`               The considered number of nearest neighbors for Uzal's
%                   L-statistic (Default is `k=3`).
% `Tw`              The considered maximum time horizon used for obtaining
%                   Uzal's L-statistic (Default is `Tw = 4*theiler`).
% `norm`            The norm used for the distance calculations. Either
%                   `euc` (Default) or `max`.
% `max_cycles`      The maximum number of embedding cycles the algorithm
%                   perform after which it breaks, no matter what (Default
%                   is `max_cycles=10`). 
%
% Output: 
%
% `Y`               The embedded phase space trajectory            
% `tau_vals`        The (different) time delays chosen by the algorithm
% `ts_vals`         The number of the chosen time series for each
%                   corresponding tau in `tau_vals`
% `epsilon_mins`    Continuity statistic. A cell array storing all epsilon-
%                   mins as a function of `delay_vals` for each encountered 
%                   time series in input `x`.
% `Ls`              The L-statistic values for the reconstruction in each
%                   encountered embedding cycle. This also includes the
%                   L-value of the very last embedding cycle, which does
%                   not contribute to the final embedding, due to an
%                   increasing L-value.

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

% define default values
delay_vals = 0:50;
sample_size = 1.0;
theiler = 1;
alpha = 0.05;
p_val = 0.5;
deltas = 13;
k = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tw = 4*theiler;  % TODO: This should be 4 times the USER-INPUT OF theiler !!!!!!!! In this case it is always 4*1=4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norm = 'euc';
expected_norms = {'euc','max'};
max_num_of_cycles = 10;

% define required and optional arguments
p = inputParser;

validScalarPosNum1 = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (x <= 1);
validScalarPosNum2 = @(x) isnumeric(x) && isscalar(x) && (x > 0) && rem(x,1)==0;
validDimension = @(x) isnumeric(x) && ismatrix(x);

addRequired(p,'x',validDimension);
addOptional(p,'delay_vals',delay_vals,validDimension);
addParameter(p,'sample_size',sample_size,validScalarPosNum1);
addParameter(p,'theiler',theiler,validScalarPosNum2);
addParameter(p,'alpha',alpha,validScalarPosNum1);
addParameter(p,'p_val',p_val,validScalarPosNum1);
addParameter(p,'deltas',deltas,validScalarPosNum2);
addParameter(p,'k',k,validScalarPosNum2);
addParameter(p,'Tw',Tw,validScalarPosNum2);
addParameter(p,'norm',norm,@(x) any(validatestring(x,expected_norms)));
addParameter(p,'max_cycles',max_num_of_cycles,validScalarPosNum2);

% parse input arguments
parse(p,x,varargin{:})

% assign variables with the resulting argument input
x = p.Results.x;
delay_vals = p.Results.delay_vals;
sample_size = p.Results.sample_size;
theiler = p.Results.theiler;
alpha = p.Results.alpha;
p_val = p.Results.p_val;
deltas = p.Results.deltas;
k = p.Results.k;
Tw = p.Results.Tw;
norm = p.Results.norm;
max_num_of_cycles = p.Results.max_cycles;

% make the input time series a column vector
if size(x,1)<size(x,2)
    x = x';
end
x_orig = x;
% normalize time series
x = (x-mean(x))./std(x);
xN = size(x,2); % number of time series

% undersampling params, which we just need for the function call of
% `pecora_embedding_cycle()`
undersamp = false;
beta = 0.05;

% Matlab in- and output check
narginchk(1,21)
nargoutchk(1,6)

%% Start computation
% set a flag, in order to tell the while loop when to stop. Each loop
% stands for encountering a new embedding dimension
flag = true;

% compute initial L-value, L_init, for all the time series and take the
% minimum
L_inits = zeros(1,xN);
for i = 1:xN
    L_inits(i) = uzal_cost(x(:,i),theiler,k,Tw,sample_size);
end
L_init = min(L_inits);

% set index-counter for the while loop
cnt = 1;

% preallocation of output
tau_vals = zeros(1,max_num_of_cycles);
ts_vals = zeros(1,max_num_of_cycles);
LS = zeros(1,max_num_of_cycles);
epsilons = cell(1,max_num_of_cycles);

% loop over increasing embedding dimensions until some break criterion will
% tell the loop to stop/break
while flag  
    % in the first run we need to find the time series to start with, i.e.
    % we have to encounter xN^2 continuity statistics...
    if cnt == 1
        % preallocate storing vectors for some variables in the coming fl
        epsilons_ = cell(1,xN);
        tau_vals__ = zeros(1,xN);
        ts_vals__ = zeros(1,xN);
        LS_ = zeros(1,xN);
        % in the multivariate case we weight each peak with the according
        % L-value
        if xN>1
            COST_ = zeros(1,xN);
        end
        
        % Loop over the potential initial time series
        for trial = 1:xN
            % call pecora_embedding_cycle for computing the continuity
            % statistic
            [epsilons_{trial}, ~, ~, ~,Y_old,~] = ...
                pecora_embedding_cycle(x,tau_vals(1:cnt),trial,...
                                   delay_vals,sample_size,theiler,...
                                   undersamp,alpha,p_val,deltas,beta,norm);

            % preallocate storing vector for the possible peaks of the 
            % continuity statistic
            tau_use = ones(1,xN); 
            Ls = zeros(1,xN);
            if xN>1
                cost = zeros(1,xN);
            end
            
            % loop over the statistics for each time series
            for ts = 1:xN
                % extract the corresponding continuity statistic
                epsilon_min_avrg = epsilons_{trial};         
                epsilon_min_avrg = epsilon_min_avrg(:,ts);
                % add 0 to account for lag 0 being also a peak
                epsilon_min_avrg = [0; epsilon_min_avrg]; 

                % Look for the local maxima:
                [pks,locs] = findpeaks(epsilon_min_avrg,'MinPeakDistance',2);

                if isempty(pks)
                    tau_use(ts)= NaN;
                    Ls(ts) = NaN;
                    continue
                else                 
                    tau_use_ = zeros(1,length(pks));
                    Ls_ = zeros(1,length(pks));
                    if xN>1
                        cost_ = zeros(1,length(pks));
                    end
                    % loop over the local maxima
                    for i = 1:length(pks)
                        % assign tau value to this peak; minus 1, because 
                        % one zero padding 
                        tau_use_(i) = delay_vals(locs(i)-1);                      

                        % compute L-statistic
                        Y_new = embed2(Y_old,x(:,ts),tau_use_(i));
                        Ls_(i) = uzal_cost(Y_new,theiler,k,Tw,sample_size);
                        if xN>1
                            cost_(i) = Ls_(i)*pks(i);
                        end
                        clear Y_new
                    end
                    
                    % pick the peak, which goes along with the least L/cost
                    if xN>1
                        [~, order] = sort(cost_);
                    else
                        [~, order] = sort(Ls_);
                    end
                    tau_use_ = tau_use_(order);
                    tau_use(ts) = tau_use_(1);
                    Ls(ts) = Ls_(order(1));
                    if xN>1
                        cost(ts) = cost_(order(1));
                        clear cost_
                    end

                    clear tau_use_, clear Ls_
                end              
            end  
              
            if xN>1               
                % Again, take the combination of time series, which yield a
                % minimum cost (=peak-height*L)
                [~,order] = sort(cost);
                tau_use = tau_use(order);             
                tau_vals__(trial) = tau_use(1);
                ts_vals__(trial) = order(1);
                LS_(trial) = Ls(order(1));
                COST_(trial) = cost(order(1));
            else
                tau_vals__(trial) = tau_use(1);
                ts_vals__(trial) = 1;
                LS_(trial) = Ls;
            end                 
        end
        
        % now compare the results for all different initial time series and
        % store the values. Sort by the smallest amount of L / cost
        if xN>1
            [~,ind] = sort(COST_);
        else
            [~,ind] = sort(LS_);
        end
        
        LS(1) = LS_(ind(1));
        ts_vals(1) = ind(1);
        ts_vals(2) = ts_vals__(ind(1)); 
        tau_vals(2) = tau_vals__(ind(1));
        epsilons{1} = epsilons_{ind(1)};
        
        % check break criteria
        if break_criteria(tau_vals, LS, max_num_of_cycles, cnt, L_init)
            tau_vals = tau_vals(1:cnt);
            ts_vals = ts_vals(1:cnt);
            LS = LS(1:cnt);
            epsilons = epsilons{1:cnt};
            flag = false;
        end
        clear Y_old
    
    % in all embedding cycles, but the first one, encounter "only" xN 
    % continuity statistics
    else
        % call pecora_embed_ts for computing the continuity statistic
        [epsilons{cnt}, ~, ~, ~, Y_old, ~] = pecora_embedding_cycle(x,...
                  tau_vals(1:cnt),ts_vals(1:cnt),delay_vals,sample_size,...
                           theiler,undersamp,alpha,p_val,deltas,beta,norm);
                
        % preallocate storing vector for the possible peaks of the continuity
        % statistic
        tau_use = ones(1,xN);
        Ls = zeros(1,xN);
        
        % loop over the statistics for each time series
        for ts = 1:xN
            % extract corresponding continuity statistic
            epsilon_min_avrg = epsilons{cnt};
            epsilon_min_avrg = epsilon_min_avrg(:,ts);
            epsilon_min_avrg = [0; epsilon_min_avrg];

            % find local maxima in the continuity statistic
            [pks,locs] = findpeaks(epsilon_min_avrg,'MinPeakDistance',2);

            if isempty(pks)
                tau_use(ts)= NaN;
                Ls(ts)=NaN;
                continue
            else
                % preallocate some variables
                tau_use_ = zeros(1,length(pks));
                Ls_ = zeros(1,length(pks));
                
                for i = 1:length(pks)
                    % now assign the tau value to this peak, minus 1, 
                    % because one zero padding
                    tau_use_(i) = delay_vals(locs(i)-1);
                    
                    % compute L-statistic
                    Y_new = embed2(Y_old,x(:,ts),tau_use_(i));
                    Ls_(i) = uzal_cost(Y_new,theiler,k,Tw,sample_size);
                    
                    clear Y_new
                end
                % pick the peak, which goes along with minimum L
                [~, order] = sort(Ls_);
                tau_use_ = tau_use_(order);

                tau_use(ts) = tau_use_(1);
                Ls(ts) = Ls_(order(1));

                clear tau_use_, clear Ls_         
            end
        end
        
        % Now pick the peak of all the "valid" ones from each time series
        if xN>1  
            % again, minimize the L-statistic
            [~,order] = sort(Ls);
            tau_use = tau_use(:,order);             
            tau_vals(cnt+1) = tau_use(1);
            ts_vals(cnt+1) = order(1);
            LS(cnt) = Ls(order(1));
        else
            tau_vals(cnt+1) = tau_use(1);
            ts_vals(cnt+1) = 1;
            LS(cnt) = Ls;
        end
        % check break criteria
        if break_criteria(tau_vals, LS, max_num_of_cycles, cnt, L_init)
            tau_vals = tau_vals(1:cnt);
            ts_vals = ts_vals(1:cnt);
            LS = LS(1:cnt);
            epsilons = epsilons{1:cnt};
            flag = false;
        end
        clear Y_old     
    end
    % increase dimension index counter  
    cnt = cnt + 1;   
end
% final reconstruction
for m = 1:length(tau_vals)
    if m == 1
        Y_tot = embed(x_orig(:,ts_vals(m)),m,tau_vals(m));
    else
        Y_tot = embed2(Y_tot,x_orig(:,ts_vals(m)),tau_vals(m));
    end
end

end

