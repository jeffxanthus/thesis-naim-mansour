% function f = ls_f(x, Prob, varargin)
%
% Computes the objective function for the least squares problem
%
%          0.5 * r(x)' * r(x)
%
% ls_f calls the TOMLAB gateway routine nlp_r to evaluate the residual r(x).

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2009 by Tomlab Optimization Inc., Sweden. $Release: 7.2.0$
% Written May 14, 1995.  Last modified Jul 17, 2009.

function f = ls_f(x, Prob, varargin)

global LS_x mad_f

if isa(x,'fmad') & Prob.ADObj == 1
    x     = getvalue(x);
    LS_x  = [];
    if nargin < 3
       r = nlp_r(x, Prob);
    else
       r = nlp_r(x, Prob, varargin{:});
    end
    f = 0.5*(r'*r);
    mad_f = f;
else
    if nargin < 3
       r = nlp_r(x, Prob);
    else
       r = nlp_r(x, Prob, varargin{:});
    end
    f = 0.5*(r'*r);
end

% MODIFICATION LOG
%
% 981023  hkh  Changed to use this routine only to compute the square sum
% 990626  hkh  Avoid feval
% 031201  hkh  Add AD handling for MAD
% 040901  med  getvalue lower case
% 060327  hkh  Use different calls to nlp_r if no varargin
% 090717  med  Fixed f calc