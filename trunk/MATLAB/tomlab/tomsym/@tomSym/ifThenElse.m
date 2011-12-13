function y=ifThenElse(cmp1,op,cmp2,yTrue,yFalse,s)
% tomSym/ifThenElse - Symbolic if/then/else
%
% y = ifThenElse(conditon, yTrue, yFalse) replicates the behaviour of
% the C-language construction y = ( condition ? yTrue : yFalse )
%
% y = ifThenElse(conditon, yTrue, yFalse, s) is a "softer" form, which
% gives a weighted average of yTrue and yFalse if the condition is close
% (within about 6*s) of being true. This makes it possible to compute the
% derivative of y with respect to variables that are used in the condition.
%
% y = ifThenElse(cmp1,op,cmp2,yTrue,yFalse,s) is a more detailed form,
% compatible with the function that is overloaded.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-01-24 by rutquist for TOMLAB release 7.7

if nargin<=4
    if islogical(cmp1)
        % The first argument had no symbols, so was evaluated to
        % true/false.
        if cmp1
            y = op;
        else
            y = cmp2;
        end
        return
    end
    if nargin<4
        yTrue = 0; % Replaces s in this construct.
    end
    if isempty(strmatch(operator(cmp1),{'eq','ne','le','lt','ge','gt'},'exact'))
        error('Condition must be a comparison.');
    end
    % Call the function with all arguments in the right place
    y = ifThenElse(operand(1,cmp1),operator(cmp1),operand(2,cmp1),op,cmp2,yTrue);
    return
elseif nargin==5
    warning('ifThenElse:wrongargs','Incorrect number of arguments to ifThenElse.')
    s = 0;
end

if size(cmp1,1)~=size(cmp2,1) || size(cmp1,2)~=size(cmp2,2)
    error('Comparison arguments are not the same size.');
end

if numel(cmp1)~=1
    siz = size(cmp1);
elseif numel(yTrue)~=1
    siz = size(yTrue);
else
    siz = size(yFalse);
end

if (numel(cmp1)~=1 && ~all(size(cmp1)==siz)) || ...
        (numel(yTrue)~=1 && ~all(size(yTrue)==siz)) || ...
        (numel(yFalse)~=1 && ~all(size(yFalse)==siz))
    error('Wrong size argument to ifThenElse');
end

if numel(cmp1)==1
    cmp1 = repmat(cmp1,siz);
end
if numel(cmp2)==1
    cmp2 = repmat(cmp2,siz);
end
if numel(yTrue)==1
    yTrue = repmat(yTrue,siz);
end
if numel(yFalse)==1
    yFalse = repmat(yFalse,siz);
end
if numel(s)==1
    s = repmat(s,siz);
end

if isequalwithequalnans(yTrue,yFalse)
    % Simplify if yTrue==yFalse
    y = yTrue;
elseif isnumeric(cmp1) && isnumeric(cmp2) && all(all(feval(op,cmp1,cmp2)))
    % Simplify if condition is always true
    y = yTrue;
elseif isnumeric(cmp1) && isnumeric(cmp2) && all(all(~feval(op,cmp1,cmp2)))
    % Simplify if condition is always false
    y = yFalse;
else
    y = tomSym(mfilename,siz(1),siz(2),cmp1,op,cmp2,yTrue,yFalse,s);
end
