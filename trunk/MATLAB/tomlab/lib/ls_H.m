% ls_H.m
%
% function H=ls_H(x, Prob, varargin)
%
% ls_H computes the Hessian to a nonlinear least squares problem
%
% First part is computed as J'*J.
% 2nd order terms are added if available
%
% ls_H calls the TOMLAB gateway routines  nlp_J and nlp_d2r

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.5.0$
% Written May 14, 1995.    Last modified Aug 14, 2006.

function H=ls_H(x, Prob, varargin)

if nargin < 3
   J = nlp_J(x, Prob);
else
   J = nlp_J(x, Prob, varargin{:});
end

H = J'*J;

if ~isempty(Prob.FUNCS.d2r)
   if nargin < 3
      H2=nlp_d2r( x, Prob, [], J);
   else
      H2=nlp_d2r( x, Prob, [], J, varargin{:});
   end
   if isempty(H2) | any(size(H2)~=size(H))
      return
   else
      H=H+H2;
   end
end

% MODIFICATION LOG
%
% 981023 hkh  Simplify this routine, only J'*J, but also new call to get
%              2nd part of Hessian, if available.
% 990626 hkh  Avoid feval
% 050121 frhe Removed r evaluation. It is done in nlp_d2r.m
% 060323 med  Call to nlp_d2r only if set
% 060327 hkh  Use different calls to nlp_J and nlp_d2r if no varargin
% 060814 med  FUNCS used for callbacks instead