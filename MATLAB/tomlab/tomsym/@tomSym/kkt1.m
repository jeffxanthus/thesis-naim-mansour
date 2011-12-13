function c=kkt1(f,c,x)
% kkt1 - First order Karush-Kuhn-Tucker conditions
%
% c = sym2prob(f,c,x) creates Karush-Kuhn-Tucker (KKT) conditons for
% minimizing f with respect to x, subject to c. (The input
% constraint set c will be a subset of the output c.)
%
% The objective f must be a tomSym symbolic object, while the constraints
% list c should be a cell array.
%
% The first order KKT conditons are necessary but not, in general,
% sufficient for optimality. The solution can be a local minumum, a local
% maximum or a saddle point. However, if f is known to be convex or have
% other special properties, then conditions given by KKT1 might be
% sufficient.
%
% KKT1 is useful in solving optimization problems where the objective
% function itself contains a convex optimization problem. This is often the
% case in game-theoretic applications. 
%
% For example, the problem:
%
%   minimize { max over y of f(x,y) subject to x>=y, y>=0 } subject to x<=2
%
% can be expressed in tomSym as
%
%   ezsolve(f(x,y),{kkt1(-f(x,y),{x>=y,y>=0},y),x<=2})
%
% assuming that -f is convex in y.
% Note that all conditons involving y are passed to KKT1 in this example.
% Feeding them directly to EZSOLVE will not give a correct solution.
%
% Inequalities in the constraints will give rise to complementary
% conditions. These are best solve using the KNITRO solver.
%
% See also: complementary

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-11-18 by rutquist for TOMLAB release 7.7

% Check indata
if numel(f)~=1
    error('The objective function must be a scalar.');
end

if ~isa(f, 'tomSym')
    error('Expected objective function to be a tomSym symbolic object.');
end

if isempty(c)
    c = {};
end
if ~iscell(c)
    c = {c};
end

% Expand list of constraints.
i = 1;
while i<=length(c)
    if iscell(c{i})
        c = {c{1:i-1}, c{i}{:} c{i+1:end}};
        continue
    end
    if ~isa(c{i}, 'tomSym')
        if ~isa(c{i},'logical')
            error(['Illegal constraint type: "' class(c{i}) ...
                '" - Constraints must be either tomSym or logical.']);
        end
        if ~all(all(c{i}))
            error('Constant constraint evaluated to "false".');
        end
        warning('tomSym:ConstConstraint','Constant constraint removed.');
        c = {c{1:i-1}, c{i+1:end}};
        continue
    end
    if strcmp(operator(c{i}),'positiveSemidefinite')
        error('Semidefinite constraints are not supported by KKT1.');
    end
    if tomCmp(c{i},'ctranspose')
        c{i} = operand(1,c{i});
    end
    if isempty(strmatch(operator(c{i}),{'eq','le','ge'}))
        error('Constraints must use one of ==, <= or >=');
    end
    if isa(operand(1,c{i}),'tomSym') && ...
            ~isempty(strmatch(operator(operand(1,c{i})),{'eq','le','ge'}))
        c = {c{1:i-1}, operand(1,c{i}), ...
            feval(operator(c{i}), operand(2,operand(1,c{i})), operand(2,c{i})), ...
            c{i+1:end}};
    elseif isa(operand(2,c{i}),'tomSym') && ...
            ~isempty(strmatch(operator(operand(2,c{i})),{'eq','le','ge'}))
        c = {c{1:i-1}, ...
            feval(operator(c{i}), operand(1,c{i}), operand(1,operand(2,c{i}))), ...
            operand(2,c{i}), c{i+1:end}};
    else
        i = i+1;
    end
end

if nargin<3
    % If x is not given, then use all symbols in f and c.
    x = {f, c{:}};
end
    
% Compile list of symbols in x (making it possible to input x in various
% ways and always end up with a vertcat).
if ~iscell(x)
    x = {x};
end
symbs = struct;
for i = 1:length(x)
    if ~isa(x{i}, 'tomSym')
        error('Expected symbol list to be a tomSym symbolic object.');
    end
    sc = symbols(c{i},'struct');
    scl = fieldnames(sc);
    for j=1:length(scl)
        symbs.(scl{j}) = sc.(scl{j});
    end
end
symbs = orderfields(symbs);
sl = fieldnames(symbs);
x = cell(size(sl));
for i=1:length(sl)
    x{i} = symbs.(sl{i});
    if operand(2,x{i})
        error('Integer variables not allowed by KKT1');
    end
    x{i} = vec(x{i});
end
x = vertcat(x{:});

% Start computing the KKT1 sum
s = derivative(f,x);
for i=1:length(c)
    switch operator(c{i})
        case {'le', 'eq'}
            g = vec(operand(2,c{i})-operand(1,c{i}));
        case 'ge'
            g = vec(operand(1,c{i})-operand(2,c{i}));
        otherwise
            error('Equations must use <=, ==, or >=.');
    end
    gx = derivative(g,x);
    if ~iszero(g)
        lam = tom([],1,numel(g));
        s = s+lam*gx;
        if ~tomCmp(c{i},'eq')
            c = {c{1:i-1}, complementary(vec(g)', lam), c{i+1:end}};
        end
    end
end
c = {c{:}, s==0};
