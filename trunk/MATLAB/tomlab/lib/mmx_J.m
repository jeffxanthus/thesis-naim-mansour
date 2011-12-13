% mmx_J.m
%
% function Jxzs=mmx_J(x, Prob, varargin)
%
% mmx_J computes the Jacobian matrix for the transformed minimax problem
% r is the residuals in the original formulation: min max r(x)
% J is the Jacobian for the original formulation
% The extra variable z=x(n+1) has derivatives -1 for the first n residuals
% and 1 for the last.
% The slack variables have the identity matrix as the Jacobian
%
% J(x,z,s)= [J,-1*e,I;0,1,0]

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Feb 5, 2000. Last modified Aug 14, 2006.

function Jxzs=mmx_J(x, Prob, varargin)

global LS_x LS_r

m = length(Prob.LS.weightY)-1;

n = length(x)-1-m;

JFunc = Prob.minimax.J;

Prob.x_0 = Prob.x_0(1:n); % Adjust x_0 for user routine
Prob.x_L = Prob.x_L(1:n);
Prob.x_U = Prob.x_U(1:n);
Prob.N   = n;

% Send also x(n+1) and the slack to p_r. 
% p_J calls user routine only with x(1:n), where n is Prob.N

if isempty(JFunc)
   rFunc=Prob.minimax.r;
   r=LS_r(1:m);
   % Always use FD algorithm
   J = FDJac(x(1:n), Prob, rFunc, r, varargin{:} ); % FD Algorithm
else
   p=xnargin(JFunc);
   if p>2
      J=feval(JFunc, x(1:n), Prob, varargin{:} );
   elseif p==2
      J=feval(JFunc, x(1:n), Prob);
   else
      J=feval(JFunc, x(1:n));
   end
end

Jxzs=sparse([J,-ones(m,1), eye(m);zeros(1,n),1,zeros(1,m)]);

% MODIFICATION LOG
%
% 000205 hkh  minimax Jacobian routine written
% 020328 hkh  modified for Tomlab v3.1
% Note - this routine is used by code infSolve, now in comments
