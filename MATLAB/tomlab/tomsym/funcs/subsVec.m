function y = subsVec(fun, x, v)
% subsVec - Expand a tomSym by substituting a vector for a scalar symbol.
%
% Y = SUBSVEC(FUN,X,V) where X is a m-by-n tomSym object, X is a scalar
% tomSym symbol, and V is a vector of length p, returns a m-by-n*p matrix
% such that Y(:,(n-1)*i+1:n*i) = subs(FUN,X,V(i)).

% TODO: Fix for subjectTo

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
    y = subsVec(subs(fun,x,lookup(v,ix)),ix,1:length(v));
    return
end

if isnumeric(fun)
    y = repmat(fun,1,length(v));
    return
end

if ~isa(fun,'tomSym')
    error(['FUN must be a tomSym. Is: ' class(fun)]);
end

if any(strcmp(operator(fun),{'eq','le','lt','gt','ge'}))
    y = feval(operator(fun), subsVec(operand(1,fun),x,v), subsVec(operand(2,fun),x,v));
    return
end

fun = expandsubjectto(fun,x,v);
if strcmp(operator(fun),'subjectTo')
   % Invert order of subsVec-subjectTo calls.
   ops = operands(fun);
   for i=3:length(ops)
       ops{i} = subsVec(ops{i},x,v);
   end
   y = subjectTo(subsVec(ops{1},x,v),ops{2},ops{3:end});
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
        y = repmat(fun,1,length(v));
        return
    end
end
    
if isempty(fieldnames(s));
    % Return numeric result
    c = cell(length(v),1);
    for i=1:length(v)
        c{i} = subs(fun,x,v(i));
    end
    y = horzcat(c{:});
else
    % Return symbolic result
    % Replace x, to avoid symbol collisions.
    ix = tom([],1,1);
    y = tomSym('subsVec',size(fun,1),size(fun,2)*length(v), subs(fun,x,ix), ix, v);   
end
