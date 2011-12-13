% nnls1  
%
% Solve non-negative least-squares by a sequence of calls to lsqnoneg
% with successively higher tolerance Tol (factor 10)
%
% function [x, w] = nnls1(A, b, Tol)
%
% INPUT:
% A     Matrix. Solve ||Ax-b||, subject to x >= 0
% b     Vector b
%
% OUTPUT:
% x     The minimizer
% w     Dual vector w where w(i) < 0 when x(i) = 0 and w(i) 
%       is approximately 0 when x(i) > 0.
% Tol   Tol is used to determine when elements of x are less than zero.
%       Default tolerance = 10*max(size(A)) * norm(A,1) * eps

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written June 22, 1999.    Last modified Oct 31, 2000.

function [x, w] = nnls1(A, b, Tol)

% Initialize Tol if not set
if nargin < 3
   Tol = 10*eps*norm(A,1)*max(size(A));
end

MaxIter=abs(log10(Tol));

optpar.TolX=Tol;
optpar.Display='off';

for i=1:MaxIter
    [x,resNorm,r,ExitFlag,Out,w]=lsqnonneg(A,b,[],optpar);
    if ExitFlag > 0, break; end
    Tol=Tol*10;
    optpar.TolX=Tol;
end