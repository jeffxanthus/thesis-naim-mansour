function idxM = getIdxM(sz)
% Get index array for a given size

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if length(sz)<2
    sz = [sz ones(1,2-length(sz))];
end

idxM = reshape(1:prod(sz),sz);
