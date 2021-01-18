function flag = break_criteria(tau_vals, LS, max_num_of_cycles, cnt, threshold)
% BREAK_CRITERIA returns a boolean, which is `true`, whenever one of the
% break criteria of the pecuzal method is met. Otherwise it returns `false`.
%
% Copyright (c) 2020
% K. Hauke Kraemer, 
% Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software and runs under MIT licence.
%%
flag = false;
% if there hasn't been any valid peak to choose from, break
if any(isnan(tau_vals))
    disp('Algorithm stopped due to non-valid time lags. NO valid embedding! (see the continuity statistic)')
    flag = true;
end

% break, if L can not be reduced anymore
if cnt == 1 && LS(cnt) > threshold
    disp('Algorithm stopped due toincreasing L-values in the first embedding cycle. NO valid embedding achieved.')
    flag = true;      
elseif cnt > 1 && LS(cnt) > threshold
    disp('Algorithm stopped due to minimum L-value reached. VALID embedding achieved.')
    flag = true;
end
% break, if max_cycles is reached
if max_num_of_cycles == cnt
    disp('Forced break of algorithm, due to reaching `max_cycles`')
    flag = true;
end   
