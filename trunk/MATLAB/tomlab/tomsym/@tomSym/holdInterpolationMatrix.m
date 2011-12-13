% holdInterpolationMatrix - Value matrix for interpPoly
% 
% Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

function M = holdInterpolationMatrix(x,xi)

if ~isnumeric(x)
    % TODO: Fix pattern/derivatives for symbolic x...
    M = tomSym(mfilename,numel(xi),numel(x),x,xi);
    return
end

x = full(x(:));
xi = vec(xi);

Mc = cell(1,length(x));

% Compute the matrix
for i = 1:length(x)
    y = zeros(size(x));
    y(i) = 1;
    pp = mkpp([x(:); Inf], y(:), 1);
    if ~isstruct(pp)
        error('TomSym feature not supported for older versions of Matlab.');
    end
    Mc{i} = ppval(pp,xi);
end

M = horzcat(Mc{:});
