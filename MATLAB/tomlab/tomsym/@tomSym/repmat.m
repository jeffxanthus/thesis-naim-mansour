function y = repmat(a,sz1,sz2,varargin)
% tomSym/repmat - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-18 by rutquist for TOMLAB release 7.7

if numel(sz1)>2 || (nargin>2 && numel(sz1)~=1)
    error('Size argument to repmat must be scalar');
end

if nargin==2;
    % Expand shorthand notation
    y = repmat(a, sz1(1), sz1(2));
else
    if numel(sz2)~=1
        error('Size argument to repmat must be scalar');
    end
    % Simplify if size = [ 1 1]
    if sz1==0 || sz2==0
        y = zeros(sz1,sz2);
    elseif sz1==1 && sz2==1
        y = a;
    else
        y = tomSym(mfilename,size(a,1)*sz1, size(a,2)*sz2,a,sz1,sz2,varargin{:});
    end
end

