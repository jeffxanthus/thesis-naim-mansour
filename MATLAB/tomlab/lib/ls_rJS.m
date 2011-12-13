% ls_rJS.m
%
% function [r_k, J_k]=ls_rJS(x, Prob, varargin)
%
% ls_rJS returns both the residual r(x) and the Jacobian J(x) for a
% nonlinear least squares problem.
% 
% ls_rJS calls the Tomlab gateway routines nlp_r, that returns the
% residual, and nlp_J, that returns the Jacobian matrix J_k.
% J_k may be sparse
%
% ls_rJS is used when calling Optimization TB 2.x.
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written July 7, 1999. Last modified Jan 29, 2003.

function [r_k, J_k]=ls_rJS(x, Prob, varargin)

nargin;
r_k = nlp_r(x, Prob, varargin{:});

if nargout > 1
   J_k = nlp_J(x, Prob, varargin{:});
end

% MODIFICATION LOG:
%
% 990707  hkh  Created from ls_rJ, but returns sparse J_k
% 030129  hkh  CURRENTLY NOT USED IN TOMLAB