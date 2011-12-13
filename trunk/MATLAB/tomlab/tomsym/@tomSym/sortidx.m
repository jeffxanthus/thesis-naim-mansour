function i = sortidx(x,dim,mode)
% tomSym/sort - Overloaded

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

if nargin < 2
    if size(x,1)==1 && numel(x)>1
        dim = 2;
    else
        dim = 1;
    end
end

if nargin < 3
    mode = 'ascend';
end

if size(x,3-dim)>1
    error('tomSym only supports sorting vectors.');
end

i = tomSym('sortidx',size(x,1),size(x,2),x,dim,mode);
