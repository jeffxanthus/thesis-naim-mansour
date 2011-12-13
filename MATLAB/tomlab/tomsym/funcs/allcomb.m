function M = allcomb(a, b, varargin)
% Allcomb - return all combinations of the rows of the input matrices
%
% M = allcomb(a,b,...)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

if nargin==1
    M = a;
elseif isempty(varargin)
    i = repmat(1:size(a,1),size(b,1),1);
    M = [a(i(:),:) repmat(b,size(a,1),1)];
else
    M = allcomb(allcomb(a,b),varargin{:});
end

