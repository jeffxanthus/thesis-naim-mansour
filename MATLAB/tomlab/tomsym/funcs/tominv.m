function x = tominv(y,x0)
% tominv - tomSym varation of inv
%
% X = TOMINV(Y) finds X so that X*Y equals the identity matrix.
%
% X = TOMINV(Y,X0) also provides a starting guess for the inverse. This may
% speed up convergence.
%
% The difference between tominv and inv is that tominv does not compute the
% inverse explicitly. Insted it creates an extra unknown and equation for
% the solver. This avoids division-by-zero and other numeric problems that
% typically arise with explicit inverses.
%
% See also: subjectTo, tommrdivide

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

if size(y,1)~=size(y,2)
    error('Matrix must be square');
end

if nargin >= 2
    x = tommrdivide(eye(size(y,2)),y,x0);
else
    x = tommrdivide(eye(size(y,2)),y);
end
