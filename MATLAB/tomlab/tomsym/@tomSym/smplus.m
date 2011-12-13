function y=smplus(a,b)
% tomSym/smplus - Scalar + matrix addition
%
% This operation maps to plus, and uses that operator, but in tomSym it
% is treated separate, because it is not the same thing as matrix + matrix
% addition.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if numel(a)~=1
    error('First argument to smplus must be a scalar.');
end

if numel(b)==1
    % Map to normal addition instead
    y = a+b;
elseif tomCmp(b,'smplus')
    % Move out scalars from b
    y = (a+operand(1,b))+operand(2,b);
elseif tomCmp(b,'vec')
    % Move addition inside of vec
    y = vec(a+operand(1,b));
elseif iszero(a)
    % Additoin by zero
    y = b;
elseif iszero(b);
    y = repmat(a,size(b));
elseif isdiag(b)
    y = setdiag(smplus(a,getdiag(b)));
else
    y = tomSym(mfilename,size(b,1),size(b,2),a,b);
end

