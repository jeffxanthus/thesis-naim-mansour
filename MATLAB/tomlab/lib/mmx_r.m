% mmx_r.m
%
% function r=mmx_r(x, Prob, varargin)
%
% mmx_r computes the minimax residual r(x,z,s) in the point x 
%
% r(x) is the residuals in the original formulation: min max r(x)
%
% r(x,z,s) = [r(x)-z*e+s;z];

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Feb 5, 2000. Last modified Apr 4, 2002.

function rxzs=mmx_r(x, Prob, varargin)

m = length(Prob.LS.weightY)-1;

n = length(x)-1-m;

rFunc = Prob.minimax.r;

Prob.x_0 = Prob.x_0(1:n); % Adjust x_0 for user routine
Prob.x_L = Prob.x_L(1:n);
Prob.x_U = Prob.x_U(1:n);
Prob.N   = n;

if xnargin(rFunc) > 2
   r = feval(rFunc, x(1:n), Prob, varargin{:} );
elseif xnargin(rFunc) > 1
   r = feval(rFunc, x(1:n), Prob );
else
   r = feval(rFunc, x(1:n));
end

rxzs=[r-x(n+1)*ones(m,1)-x(n+2:n+m+1);x(n+1)];

% MODIFICATION LOG
%
% 000205 hkh  minimax residual routine written
% 020328 hkh  modified for Tomlab v3.1
% Note - this routine is used by code infSolve, now in comments
