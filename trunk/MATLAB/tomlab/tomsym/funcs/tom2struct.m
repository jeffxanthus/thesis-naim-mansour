function g = tom2struct(g,eq)
% tom2struct - Convert tomSym equations into a struct
%
% g = tom2struct(eq) converts eq to a struct.
% g = tom2struct(g,eq) adds the tomSym equation eq to the structure g,
%     substituting any expressions already in the struct into the result.
% g = tom2struct(g,{eq1, eq2}) converts, in turn, eq1 and eq2.
%
% The equation EQ should have the form of a tomSym equality, where one side
% is a simple symbol, and the other side is a tomSym expression or a
% constant. Alternatively, EQ may be a linear equation involving a single
% new variable and any number of previously defined variables.
%
% TOM2STRUCT does variable substitutions for any variables that are already
% defined in G, but it does not solve systems of equations. It is therefore
% important to provide equations in the correct order, so that each new
% equation only contains one new variable.
%
% Example:
% tom2struct({x == 2, 2*y == x}) is the same as struct('x',2,'y',1)
%
% tom2struct is used by ezsolve to process the starting guess.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-08-02 by rutquist for TOMLAB release 7.7

if nargin==1
    eq = g;
    g = struct;
end

if iscell(eq)
    eq = flattencellarray(eq);
    shrt = false(1,length(eq));
    for i=1:numel(eq)
        if tomCmp(eq{i},'eq') && length(symbols(eq{i})) == 1
            g = tom2struct(g,eq{i});
            shrt(i) = true;
        end
    end
    for i=1:numel(eq)
        if ~shrt(i)
            g = tom2struct(g,subs(eq{i},g));
        end
    end
    return
end

if isstruct(eq)
   fn = fieldnames(eq) ;
   for i=1:length(fn)
       g.(fn{i}) = eq.(fn{i});
   end
   return
end

if ~tomCmp(eq,'eq')
    disp(eq)
    error('Expected a tomSym equation');
end

eq_orig = eq;
eq = subs(eq,g);

if ~isa(eq,'tomSym')
    disp('Check that your list does not contain more equations than unknowns.');
    error(['Equation did not contain any new symbol: ' mcodestr(eq_orig)]);
end

lhs = operand(1,eq);
rhs = operand(2,eq);

if tomCmp(lhs,'tom') && ~isdependent(rhs,lhs)
    if tomCmp(rhs,'srepmat')
        rhs = double(rhs);
    end
    if isa(rhs,'tomSym')
        warning('tom2struct:notConstant',...[
            ['The right-hand-side of equation is not a constant: ' mcodestr(eq)]);
    end
    g.(operand(1,lhs)) = rhs;
elseif tomCmp(rhs,'tom') && ~isdependent(lhs,rhs)
    if isa(lhs,'tomSym')
        warning('tom2struct:notConstant',...
            ['The left-hand-side of equation is not a constant: ' mcodestr(eq)]);
    end
    g.(operand(1,rhs)) = lhs;
else
    f = rhs-lhs;
    x = symbols(f,'vector');
    J = derivative(f,x);
    if ~isnumeric(J)
        disp('You may only use linear equations to define a starting guess.');
        error(['tom2struct does not solve nonlinear equations: ' mcodestr(eq)]);
    end

    %TODO: Error/warning when solution seems wrong.
    %if rank(full(J)) ~= length(x)
    %    error(['Equation does not have a unique solution:' mcodestr(eq)]);
    %end
    
    f0 = constpart(f);
    x0 = -J\vec(f0);
    gg = symbols(x,'struct');
    j = 0;
    if tomCmp(x,'vertcat')
        for i = 1:nOperands(x)
            s = symbols(operand(i,x));
            g.(s{1}) = reshape(x0(j+1:j+numel(gg.(s{1}))),size(gg.(s{1})));
            j = j+numel(gg.(s{1}));
        end
    else
        if tomCmp(x,'vec')
            x = operand(1,x);
        end
        g.(operand(1,x)) = reshape(x0,size(x));
    end
end
