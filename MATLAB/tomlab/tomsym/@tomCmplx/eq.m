function o = eq(a,b)
% tomCmplx/eq - Comparison of complex-valued object
%
% Returns a cell array of two tomSym equations, unless one of these
% evaluates to a true or false, in which case a logical or tomSym is
% returned.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

o = {real(a) == real(b), imag(a) == imag(b)};

% Simplify if one of the equations evaluated to true/false.
if islogical(o{1})
    if ~o{1}
        o = o{1};
    else
        o = o{2};
    end
elseif islogical(o{2})
    if ~o{2}
        o = o{2};
    else
        o = o{1};
    end
end

