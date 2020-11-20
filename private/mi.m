function varargout=mi(varargin)
%MI   Histogram based mutual information.
%    I=MI(X) computes the auto mutual information using 10 bins,
%    where X can be a multi-column vector.
%
%    I=MI(X1,X2,...,Xn) computes the pair-wise mutual information
%    between all pairs of X1,X2,...,Xn using 10 bins. The result I
%    is a N x N matrix.
%
%    I=MI(X1,X2,...,Xn,L) or I=MI(X1,X2,...,Xn,N,L), where L is a 
%    scalar, computes mutual informations for shifting of the 
%    components of each pair of X1,X2,...,Xn until a maximal lag L. 
%
%    I=MI(X1,X2,...,Xn,N,L), where N is a scalar determining the number
%    of bins N.
%
%    [I S]=MI(...) computes the mutual information and the 
%    standard error (only for one- and two-dimensional data).
%
%    Remark
%    Please note that the mutual information derived with MI slightly 
%    differs from the results derived with MIGRAM. The reason is that
%    MI also considers estimation errors. 
%
%    Examples: x = sin(0:.2:8*pi)' + .1*randn(126,1);
%              mi(x,40,100)
%
%    References: 
%    Roulston, M. S.:
%    Estimating the errors on measured entropy and mutual 
%    information, Physica D, 125, 1999.

% Copyright (c) 2008-2020
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the input

narginchk(1,99)
nargoutchk(0,2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read the input

lag=[]; 
nbin=10;
x={}; sm_length=[];
if nargin && isnumeric(varargin{1})

  sum_m=0; 

  for i=1:nargin
    if isnumeric(varargin{i}) && max(size(varargin{i}))>1
      y=varargin{i};
      if size(y,1)==1, y=y'; end 
      my=size(y,2); sum_m=sum_m+my;
      mx=size(y,1); 
      if isempty(sm_length), sm_length=mx; end
      if mx < sm_length, sm_length=mx; end
      x(i)={y};
    end
  end
  m=size(x,2); if sum_m==1; lag=0; x(2)=x(1); m=2; end
  if isempty(x), error('Not a valid input vector.'), end

  i_double=find(cellfun('isclass',varargin,'double'));
  i_char=find(cellfun('isclass',varargin,'char'));

  if length(i_double)>1
    if max(size(varargin{i_double(end-1)}))==1
      nbin=varargin{i_double(end-1)}; 
    end
  end
  if max(size(varargin{i_double(end)}))==1
    lag=varargin{i_double(end)}; else if isempty(lag), lag=0; end
  end

  	      
else
  error('Input required')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% computation 


  MI = zeros(m,m,lag);
  MI_sigma = zeros(m,m,lag);

  for t=0:lag,
  
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the distributions

     for k=1:m; 
         for l=1:m;
             P2=hist2(x{k},x{l},nbin,t,'silent');    % 2D distribution
             mP=length(size(P2));
             %disp(['k: ',num2str(k),' l: ',num2str(l),' test: ',num2str(sum(~P2(:))/numel(P2))])
             if sum(~P2(:))/numel(P2) > 0.9995
                 warning on  
                 warning('Too less data points for the estimation of the joint distribution.'); 
                 MI(k,l,(t:lag)+1)=NaN;
                 MI_sigma(k,l,(t:lag)+1)=NaN;
                 %MI(k,l,(t)+1)=NaN;
                 %MI_sigma(k,l,(t)+1)=NaN;
                 if nogui~=2, delete(hw); end
                 return
             end
             P2=P2/sum(P2(:));      % normalization

             clear P1x
             if mP==2
                 P1x=sum(P2,2);           % distribution for x
                 P1y=sum(P2);             % distribution for y
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the mutual information
                 I1=[-sum((P1x(P1x~=0)).*log(P1x(P1x~=0))), -sum((P1y(P1y~=0)).*log(P1y(P1y~=0)))];                      % entropies of Px and Py
                 I2=-sum(P2(P2~=0).*log(P2(P2~=0)));                                                                    % entropy of joint distribution Pxy
                 I2_syserr=(length(P1x(P1x~=0))+length(P1y(P1y~=0))-length(P2(P2~=0).*log(P2(P2~=0)))-1)/(2*length(x{k})); % standard error for estimation of entropy
                 MI(k,l,t+1)=I1(1)+I1(2)-I2+I2_syserr;
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the standard errors
                 P2xy=P1x*P1y;
                 i=(P2xy~=0 & P2~=0);
                 MI_sigma(k,l,t+1)=sqrt( sum( (log(P2xy(i)./P2(i))+MI(t+1)).^2 .*(P2(i).*(1-P2(i)))) /length(x{k}));
             else
                 I1=[];
                 for i=1:mP;
                     P21=permute(P2,[i,i+1:mP, 1:(i-1)]);
                     P1x(:,i)=sum(reshape(P21,size(P21,2),size(P21,2)^(mP-1)),2);   % distribution for x_i
                     I1 = -sum((P1x(P1x~=0)).*log(P1x(P1x~=0))); % sum of entropies of Px_i
                 end  
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the mutual information
                 I2=-sum(P2(P2~=0).*log(P2(P2~=0)));                                                                    % entropy of joint distribution Pxy
                 I2_syserr=(length(P1x(P1x~=0))-length(P2(P2~=0).*log(P2(P2~=0)))-1)/(2*length(x{k})); % standard error for estimation of entropy
                 MI(k,l,t+1)=I1-I2+I2_syserr;
             end
         end
     end
     
         
  end


  if nargout==1
     varargout(1)={MI};
  elseif nargout==2
     varargout(1)={MI};
     varargout(2)={MI_sigma};
  end

  warning on
  if isnan(MI(1,1)), warning('Mutual information is NaN. Use smaller bin size.'), end
  try set(0,props.root), catch end

  if nargout==0 && nogui
      varargout{1} = MI;
  end


        

