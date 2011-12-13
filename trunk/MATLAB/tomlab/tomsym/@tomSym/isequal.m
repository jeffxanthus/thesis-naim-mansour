function y=isequal(a,b)
% tomSym/isequal - Overloaded function
%
% In tomSym/isequal different ids do not result in "false" if the objects
% are otherwise identical.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if ~strcmp(class(a),class(b))
    y = false;
    return
end

if isequal(a.s(end).id,b.s(end).id)
    y = true;
    return
end

if ~isequal(a.d, b.d) || length(a.s) ~= length(b.s)
    y = false;
    return
end


for i=1:length(a.s)
    if ~(isequal(a.s(i).op, b.s(i).op) && isequal(a.s(i).a, b.s(i).a))
        y = false;
        return
    end
end

y = true;
