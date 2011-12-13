function y = isint(x)
% tomSym/isint - Determine if a tomSym is guaranteed integer
%
% y = isint(x) returns "true" if x contains only integer symbols and
% constants, and uses only operations that are guranteed to yield integer
% results.
%
% If "false" is returned, this is no guarantee that x is non-integer. It
% just means that isinteger was unable to prove it to be integer.
%
% isint should not be confused with Matlab's builtin "isinteger". The
% latter will return "false" for any tomsym.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7


y = false;

for i=1:length(x.d)
    if isnumeric(x.d{i}) && ~isint(x.d{i})
        return
    end
end

for i=1:length(x.s)
    switch x.s(i).op
        case 'tom'
            if ~x.d{-x.s(i).a(2)}
                return
            end
        case {'plus','minus','uplus','uminus','times','mtimes','power',...
                'mpower','spower','transpose','ctranspose','horzcat',...
                'vertcat','sparse','setSymmetric','reshape','repmat',...
                'vec','lookup','submatrix'}
            % Do nothing
        otherwise
            return
    end
end

y = true;
