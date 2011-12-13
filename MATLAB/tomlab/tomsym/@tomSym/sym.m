function s = sym(o)
% tomSym/sym - convert a tomSym object to a symbolic toolbox object
%
% s = sym(o) creates a symbolic object s from the tomSym object o.
%
% This makes it possible to use symbolic toolbox features to manipulate
% a tomSym.
%
% If all the symbols used in the expression are tomSym symbols in the
% current workspace (i.e. if they were declared using toms) then the
% inverse transformation can be accomplished by o = eval(s).

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

arg = o.arguments;
for i=1:length(arg)
    if isa(arg{i}, 'tomSym')
        arg{i} = sym(arg{i});
    end
end

switch(operator(o))
    case 'constant'
        s = arg{1};
    case 'tom'
        s = sym(arg{1});
    case 'function'
        s = feval(arg{:});
    otherwise
        s = feval(operator(o), arg{:});
end
