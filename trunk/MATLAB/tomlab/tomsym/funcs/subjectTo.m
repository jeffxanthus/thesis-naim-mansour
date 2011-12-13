function y = subjectTo(f,x0,varargin)
% subjectTo - Create a combined function/constraint object for tomSym
%
% y = subjectTo(f,guess,c1,c2,...) returns an object that evaluates to f, but
% carries with it the constraints c1, c2, .... This is useful, for example, 
% when a function returning a tomSym object needs to pass contraints along 
% with it.
%
% X0 should be a struct or cell array listing starting guesses for the 
% free variables. The remaining variables are considered parameters. 
% If there are no parameters in the constraints (i.e. x lists every
% variable that is used in c1, c2, ...) then subjectTo will solve for the
% unknowns. Otherwise it is treated as a symbolic expression.
%
% The free variables will be unique when substitutions are made 
% (for example when using collocation with PROPT.)
%
% Example: 
% The following is a way of avoiding the sharp nonlinearity in y=abs(x)
% which works if the objective involves minimizing y.
% s = tom([],size(x,1),size(x,2)) % anonymous symbol 
% y = subjectTo(s,s==99,-s <= x, x <= s)
%
% See also: tomSym/fzero

% (See also: extractConstraints)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

c = varargin;
s = symbols(tomSym(c),'struct');

if iscell(x0) || isa(x0,'tomSym')
    ws = warning('off','tom2struct:notConstant');
    x0 = tom2struct(x0);
    warning(ws);
end

p = s;
fn = fieldnames(x0);
for i=1:length(fn)
    if ~isfield(p,fn{i})
        warning('subjectTo:unusedx0', ...
            ['Initial guess for ' fn{i} ' which is not used.']);
        x0 = rmfield(x0, fn{i});
    else
        p = rmfield(p,char(fn{i}));
    end
end

if isempty(fieldnames(p))
    % No parameters -  we can actually compute y
    options = struct;
    options.prilev = 0;
    sol = ezsolve(0,c,x0,options);
    y = subs(f,sol);
else
    y = tomSym(mfilename, size(f,1), size(f,2), f, x0, varargin{:});
end
