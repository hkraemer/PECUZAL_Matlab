function [distances, comp_dist] = all_distances(varargin)
% ALL_DISTANCES Distance between given point to other points in phase space.
%    D = ALL_DISTANCES(FID_POINT,Y) computes the vector D of length N
%    of distances between the M-size vector FID_POINT and the M-size 
%    vectors given in NxM-matrix Y.
% 
%    [D, C] = ALL_DISTANCES(FID_POINT,Y) computes the matrix C of size
%    NxM of distances between point FID_POINT and the points in Y for 
%    each of the M components separately.
%
%    ... = ALL_DISTANCES(FID_POINT,Y,NORM) specifies the norm used to
%    calculate distances:
%
%      'euc' - (default) Euclidean norm
%      'max' - maximum norm
%
%    This function is meant to determine the neighbourhood of a certain point
%    without computing the whole distances matrix (as being done by the 
%    PDIST()-function).

% Copyright (c) 2020
% K. Hauke Kraemer, 
% Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software and runs under MIT licence.

fid_point = varargin{1};
Y = varargin{2};

methLib={'euc', 'max'}; % the possible norms
try
    meth = varargin{3};
    if ~isa(meth,'char') || ~ismember(meth,methLib)
       warning(sprintf('Specified norm should be one of the following possible values: ''%s'', ''%s''.', methLib{:}))
    end
catch
    meth = 'euc';
end

YY = repmat(fid_point, size(Y,1), 1);
if strcmp(meth, 'euc')
    distances = sqrt(sum((YY - Y) .^ 2, 2));
elseif strcmp(meth, 'max')
    distances = max(abs(YY - Y), [], 2);
end

if nargout > 1
    comp_dist = abs(YY - Y);
end
