function Y = scalerows(sc,M)
% tomSym/scalerows - Overloaded function
%

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-09-08 by rutquist for TOMLAB release 7.7

if size(sc,1) ~= size(M,1) || size(sc,2)~=1
    error('Scale vector must be a column vector with the same number of rows as M');
end

% Move out minus signs
if tomCmp(sc,'uminus')
    Y = -scalerows(-sc,M);
    return
end
if tomCmp(M,'uminus')
    Y = -scalerows(sc,-M);
    return
end

% Change repmat into smtimes
if tomCmp(sc,'repmat') && numel(operand(1,sc))==1
    Y = smtimes(operand(1,sc),M);
    return
end

if size(M,1)==1
    Y = sc.*M; % Scalar times row
elseif size(M,2)==1
    Y = sc.*M; % element-wise product
elseif isnumeric(sc) && all(diff(sc)==0)
    Y = sc(1).*M;
elseif tomCmp(M,'scalerows')
    Y = scalerows(sc.*operand(1,M),operand(2,M));
elseif isdiag(M)
    Y = setdiag(sc.*getdiag(M));
else
    Y = tomSym(mfilename,size(M,1),size(M,2),sc,M);
end
