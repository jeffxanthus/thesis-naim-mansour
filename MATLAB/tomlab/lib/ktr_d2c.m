% TOMLAB gateway routine
% Callback from KNITRO
%
% ktr_d2c computes the Hessian
%
%       lam * d2c(x)
%
% of the nonlinear constraints
%
% function d2c = ktr_d2c(x, lam, Prob)

% Anders Göran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Nov 25, 2011.   Last modified Nov 25, 2011.

function d2c = ktr_d2c(x, lam, Prob)

global n_H n_d2c NARG

x=x(:);

if(isempty(lam))
   % If no nonlinear c/s, should not really happen
   z = spalloc(Prob.N,Prob.N,0);
else
   z = nlp_d2c(x,lam,Prob);
end

d2c = triu(z);
if ~issparse(d2c), d2c = sparse(d2c); end

% MODIFICATION LOG:
%
% 111125 ango Wrote file, based on ktr_d2L.m
