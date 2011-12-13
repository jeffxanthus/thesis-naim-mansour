function v=getdiag(M,d)
% tomSym/getdiag - Getdiag for tomSym.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if size(M,1)~=size(M,2) || nargin>1
    % Get off-center diagonals
end

% Simplify if the argument is zero or scalar.
if tomCmp(M,'constant')
    v = getdiag(operand(1,M));
elseif numel(M)==1
    v = M;
elseif tomCmp(M,'setdiag')
    v = vec(operand(1,M));
else
    v = tomSym(mfilename,size(M,1),1,M);
end
