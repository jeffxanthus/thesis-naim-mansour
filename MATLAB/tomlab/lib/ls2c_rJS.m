% ls2c_rJS.m
%
% function [r_k, J_k]=ls2c_rJS(x, Prob, varargin)
%
% ls2c_rJS is used for Curve Fitting problems
%
% ls2c_rJS returns both the residual r(x) and the Jacobian J(x) for a
% nonlinear least squares problem.
% J_k may be sparse
%
% ls2c_rJS calls the user function that returns the
% residual, and possibly the user function that returns the Jacobian matrix.
%
% ls2c_rJS is used when implementing the OPTIM TB 2.0 compatibility interface

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written July 28, 1999. Last modified Aug 14, 2006.

function [r_k, J_k]=ls2c_rJS(x, Prob, varargin)

nargin;
global n_r n_J

n_r = n_r + 1;
r=Prob.FUNCSX.r;

if xnargin(r) > 2
   r_k = feval(r, x, Prob);
else
   r_k = feval(r, x, Prob.LS.t);
end
r_k=r_k-Prob.LS.y;

if nargout > 1
   J=Prob.FUNCSX.J;
   if ~isempty(J)
      n_J = n_J + 1;
      if xnargin(J) > 2
         J_k = feval(J, x, Prob);
         %J_k = feval(J, x, Prob.LS.t, Prob);
      else
         J_k = feval(J, x, Prob.LS.t, Prob);
      end
   else
      J_k=[];
   end
end

% MODIFICATION LOG:
%
% 060814  med  FUNCS used for callbacks instead