function y = fsum(fun, x, v)
% fsum - Expand a tomSym by substituting a vector for a scalar symbol.
%
% Y = FSUM(FUN,X,V) where X is a tomSym object, X is a scalar tomSym
% symbol, and V is a vector, returns the sum of FUN for X taking the values
% in V. Besides efficiency, FSUM equivalent to the code code:
%
% Y = 0
% for i=1:length(V)
%    Y = Y+subs(FUN,X,V(i))
% end
%
% FSUM is more efficient than the equivalent for loop, because the for-loop
% creates a new symbolic term for each element of V, which can result in a
% very long expression for Y.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-10 by rutquist for TOMLAB release 7.7

if ~(tomCmp(x,'tom') && numel(x)==1)
    error('X must be a scalar symbol.');
end

if length(v) ~= numel(v)
    error('V must be a vector.');
end

if ~isnumeric(v) || v(end)~=length(v) || ~all(diff(v)==1)
    % Re-write so that v is always 1:N
    ix = tom([],1,1);
    y = fsum(subs(fun,x,lookup(v,ix)),ix,1:length(v));
    return
end

fun = expandsubjectto(fun,x,v);
if strcmp(operator(fun),'subjectTo')
   % Invert order of fsum-subjectTo calls.
   ops = operands(fun);
   for i=3:length(ops)
       ops{i} = subsVec(ops{i},x,v);
   end
   y = subjectTo(fsum(ops{1},x,v),ops{2},ops{3:end});
   return
end

s  = symbols(fun,'struct');
if isfield(s,operand(1,x))
    s = rmfield(s,operand(1,x));
else
    if isfield(symbols(fun,'all'),operand(1,x))
        error(['When using subsVec/fsum inside another subsVec/fsum, ' ...
            'the index names must be different.']);
    else
        y = length(v)*fun;
        return
    end
end
    
if isempty(fieldnames(s));
    % Return numeric result
    y = zeros(size(fun));
    for i=1:length(v)
        y = y+subs(fun,x,v(i));
    end
else
    % Return symbolic result
    % Replace x, to avoid symbol collisions.
    ix = tom([],1,1);
    y = tomSym('fsum',size(fun,1),size(fun,2), subs(fun,x,ix), ix, v);   
end
