% function [F, J, domerr] = mcp_funjac(x, jacflag, Prob)
%
% TOMLAB gateway routine for the computation of PATH information.
%
% INPUT:
% x       Current iterate x.
% jacflag Indicates whether or not the Jacobian is needed.
% Prob    TOMLAB problem structure.
%
% To compute F(x) nlp_r calls the routine Prob.FUNCS.r as 
%            r=feval(Prob.FUNCS.r, x)
%
% The global counter variable n_r is incremented
%
% To compute J(x) nlp_J calls the routine Prob.FUNCS.J as 
%            J=feval(Prob.FUNCS.J, x)
% 
% OUTPUT:
% F       Function evaluations at x.
% J       The Jacobian of F.
% domerr  The number of divided by zero (set to 0).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Apr 18, 2004.  Last modified Apr 18, 2004.

function [F, J, domerr] = mcp_funjac(x, jacflag, Prob)

F = [];
J = [];
domerr = 0;

F = nlp_r(x,Prob);
if jacflag
    J = sparse(nlp_J(x,Prob));
end