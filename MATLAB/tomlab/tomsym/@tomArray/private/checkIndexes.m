function checkIndexes(o)
% Check that index names of tomArray are OK.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

named = false(size(o.sz));
for i=1:length(o.ni)
    named(i) = ~isempty(o.ni{i});
end
if ~(isempty(o.ni) || all(named))
    error('Either all indexes must be named, or none.');
end

if length(unique(o.ni))~=length(o.ni)
    error('Index names must be unique');
end
