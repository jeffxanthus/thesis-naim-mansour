function Y = scalecolumns(sc,M)
% tomSym/scalecolumns - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-08-24 by rutquist for TOMLAB release 7.7

if length(sc) ~= size(M,2) || size(sc,2)~=1
    error('Scale must be a column vector of length equal to the number of columns of M');
end

% Move out minus signs
if tomCmp(sc,'uminus')
    Y = -scalecolumns(-sc,M);
    return
end
if tomCmp(M,'uminus')
    Y = -scalecolumns(sc,-M);
    return
end

if size(M,1)==1
    Y = sc'.*M; % element-wise product
elseif size(M,2)==1
    Y = sc.*M; % scalar times column
elseif tomCmp(M,'scalecolumns')
    Y = scalecolumns(sc.*operand(1,M),operand(2,M));
elseif isdiag(M)
    Y = setdiag(sc.*getdiag(M));
else
    Y = tomSym(mfilename,size(M,1),size(M,2),sc,M);
end
