% function g=ls_g(x, Prob, varargin)
%
% Computes the gradient to a nonlinear least squares problem,
%
%            J(x)' * r(x)
%
% ls_g calls the TOMLAB gateway routines nlp_r, that returns the
% residual r(x), and nlp_J, that returns the Jacobian matrix J(x).

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2009 by Tomlab Optimization Inc., Sweden. $Release: 7.3.0$
% Written May 14, 1995.  Last modified Aug 13, 2009.

function g = ls_g(x, Prob, varargin)

global mad_J mad_r LS_x mad_g

if isa(x,'fmad') & Prob.ADObj == -1
    x = getvalue(x);
    LS_x = []; LS_xJ = [];
    Prob.ADObj = 1;
    if nargin < 3
       nlp_r(x, Prob);
       Prob.ADObj = -1;
       nlp_J(x, Prob);
    else
       nlp_r(x, Prob, varargin{:});
       Prob.ADObj = -1;
       nlp_J(x, Prob, varargin{:});
    end
    g     = mad_J'*mad_r;
    mad_g = g;
else
    if nargin < 3
       r = nlp_r(x, Prob); 
       J = nlp_J(x, Prob);
    else
       r = nlp_r(x, Prob, varargin{:}); 
       J = nlp_J(x, Prob, varargin{:});
    end
    if isempty(r) || isempty(J)
        g = [];
    else
        g = J' * r;
    end
end

% MODIFICATION LOG
%
% 981023  hkh  Simplify this routine, only J'*g. Call nlp_r/nlp_J instead
% 990626  hkh  Avoid feval
% 031201  hkh  Adding AD handling for MAD
% 040901  med  getvalue lower case
% 060327  hkh  Use different call to nlp_r and nlp_J if no varargin
% 090813  med  mlint check