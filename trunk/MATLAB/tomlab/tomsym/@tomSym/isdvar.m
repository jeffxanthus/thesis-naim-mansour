function y = isdvar(x)
% tomSym/isdvar - Determine if a tomSym is a decision variable
%
% y = isdvar(x) returns "true" if x contains only symbols that translate
% directly to decision variables. That is: no mathematical expressions are
% allowed. (The only allowed operations are things like transpose, index
% lookup and reshape.)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

y = false;

for i=1:length(x.s)
    switch x.s(i).op
        case {'tom','transpose','ctranspose','horzcat','lookup','submatrix'...
                'vertcat','sparse','setSymmetric','reshape','repmat','vec'}
            % Do nothing
        otherwise
            return
    end
end

y = true;
