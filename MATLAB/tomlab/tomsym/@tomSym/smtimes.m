function y=smtimes(a,b)
% tomSym/smtimes - Scalar x matrix multiplication
%
% This operation maps to mtimes, and uses that operator, but in tomSym it
% has to be treated separate, because it is not the same thing as matrix
% multiplication. (It is the same thing as kron, but is kept separate from
% that too, in order to be displayed as in matlab.) The scalar is always
% moved to the left.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-06-15 by rutquist for TOMLAB release 7.7

if numel(a)~=1
    error('First argument to smtimes must be a scalar.');
end

% Check for complex constants
if isnumeric(a) && ~isreal(a) && any(imag(a(:)))
   warning('tomsym:complexdata','Complex data encountered. Attempting to convert to tomCmplx.');
   b = tomCmplx(b,spzeros(size(b)));
   y = feval('mtimes',a,b);
   return
end
if isnumeric(b) && ~isreal(b) && any(imag(b(:)))
   warning('tomsym:complexdata','Complex data encountered. Attempting to convert to tomCmplx.');
   a = tomCmplx(a,spzeros(size(a)));
   y = feval('mtimes',a,b);
   return
end

if numel(b)==1
    % Map to normal multiplication instead
    y = a*b;
elseif tomCmp(b,'smtimes')
    % Move out scalars from b
    y = (a*operand(1,b)).*operand(2,b);
elseif tomCmp(b,'vec')
    % Move multiplication inside of vec
    y = vec(a*operand(1,b));
elseif tomCmp(b,'setdiag')
    % TODO: There are more functions that should be treated this way.
    y = setdiag(a*operand(1,b));
elseif iszero(a) || iszero(b)
    % Multiplication by zero
    y = spzeros(size(b));
elseif isone(a)
    % Multiplication by 1
    y = b;
elseif isnumeric(a) && a==-1
    % Multiplication by -1
    y = -b;
elseif isallone(b)
    y = repmat(a,size(b));
elseif isdiag(b)
    y = setdiag(smtimes(a,getdiag(b)));
else
    y = tomSym(mfilename,size(b,1),size(b,2),a,b);
end
