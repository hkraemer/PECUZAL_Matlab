function varargout=hist2(varargin)
%HIST2   Two-dimensional histogram.
%    P=HIST2(X) bins the two-dimensional density of X(i) and 
%    X(i+1) into a 10x10 equally spaced matrix and returns it 
%    in P.
%
%    P=HIST2(X,Y) bins the two-dimensional density of X(i) and 
%    Y(i) into a 10x10 equally spaced matrix and returns it 
%    in P.
%
%    P=HIST2(X,K,L), where K and L are scalars, uses K bins and
%    a lag L.
%
%    [P,J]=HIST2(...) returns the matrix P and the two-dimensional
%    vector J containing the two-dimensional density matrix and the
%    bin location for X (and Y).
% 
%    HIST2(...) without any output arguments produces a histogram plot.
%
%    Examples: x = randn(10000,1);
%              hist2(x)
%
%    See also HIST, HISTN, BAR3, MI.

% Copyright (c) 2008-2020
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or any later version.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check the input

narginchk(1,5)
nargoutchk(0,2)
varargin{6}=[];
lag=[]; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read the input

x=varargin{1};
if isempty(varargin{2}) | ~isnumeric(varargin{2}), y=x; lag=[];
else, y=double(varargin{2}); end

i_double=find(cellfun('isclass',varargin,'double'));
i_char=find(cellfun('isclass',varargin,'char'));

if (isnumeric(varargin{2}) & max(size(varargin{2}))==1) | ~isnumeric(varargin{2})
   y=x;
   if ~isempty(varargin{i_double(2)}), nbin=varargin{i_double(2)}; else nbin=10; end
   if ~isempty(varargin{i_double(3)}), lag=varargin{i_double(3)}; else lag=1; end
else
   if ~isempty(varargin{i_double(3)}), nbin=varargin{i_double(3)}; else nbin=10; end
   if ~isempty(varargin{i_double(4)}), lag=varargin{i_double(4)}; else if isempty(lag), lag=0; end, end
end


if ~(prod(double(size(y)==size(x)))), error('Both vectors must have the same size.'), end
if size(x,1)<size(x,2), x=x';end
if size(y,1)<size(y,2), y=y';end
x=x(:,1); y=y(:,1);
%x=x(:); y=y(:);

% normalise the value range to [0 1)
x1 = (x - min(x)) / max(x - min(x)) - eps;
y1 = (y - min(y)) / max(y - min(y)) - eps;

if length(x)<=lag, lag=length(x)-1; warning(['Lag is too large. Using ',num2str(lag),' instead.']), end

% this is the main trick: put the first data to values [1:1:nbin] and the second data to [nbin:nbin:nbin^2]
temp = fix(x1(1:(end-lag)) * nbin) + (fix(y1((lag+1):end) * nbin)) * nbin;

% call Matlab histc function (faster than hist)
p=histc(temp,0:(nbin^2-1))';
p2=reshape(p,nbin,nbin)/length(x);

% create the correct bin denominator
minx = min(min(x)); maxx = max(max(x));
if minx == maxx, minx = minx - floor(nbin/2) - 0.5; maxx = maxx + ceil(nbin/2) - 0.5; end
miny = min(min(y)); maxy = max(max(y));
if miny == maxy, miny = miny - floor(nbin/2) - 0.5; maxy = maxy + ceil(nbin/2) - 0.5; end

bins(1,:)=minx:(maxx-minx)/(nbin-1):maxx+eps;
bins(2,:)=miny:(maxy-miny)/(nbin-1):maxy+eps;
if nargout==0
  p2 = flipud(p2);
  bar3(bins(2,:),p2')
  dhy=bins(2,2)-bins(2,1);
  set(gca,'xlim',[0 nbin+1-.1],'ylim',[bins(2,1)-dhy bins(2,end)+dhy])
  hx=get(gca,'xtick'); dhx=diff(hx); 
  set(gca,'xtick',[1:dhx(1):nbin+1])
  if dhx>1
    tx=num2str(fix(100*[bins(1,1):(bins(1,end)-bins(1,1))/(length(hx)-1):bins(1,end)]')/100);
  else
    tx=num2str(fix(100*[bins(1,:)]')/100);
  end
  set(gca,'xticklabel',tx)
  xlabel('x'), ylabel('y')
elseif nargout==1
  varargout{1}=fliplr(p2);
elseif nargout==2
  varargout{1}=fliplr(p2);
  varargout{2}=bins';
end
