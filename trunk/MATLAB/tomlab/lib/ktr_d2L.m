% TOMLAB gateway routine
% Callback from KNITRO
%
% ktr_d2L computes the Hessian
%
%       d2f(x) + lam * d2c(x)
%
% to the Lagrangian function,
%
%   L(x,lam) =   f(x) + lam' * c(x)
%
% function d2L=ktr_d2L(x, lam, Prob)

% Anders Göran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written July 8, 2003.   Last modified Aug 13, 2009.

function d2L=ktr_d2L(x, lam, Prob)

global n_H n_d2c NARG

x=x(:);

if(isempty(lam))
   % If no nonlinear c/s
   z = nlp_H(x,Prob);
else
   z = nlp_H(x,Prob);
   if isempty(z)
      z = nlp_d2c(x,lam,Prob);
   else
      % Fix for AD to work
      v=nlp_d2c(x,lam,Prob);
      if size(v,1) ~= size(z,1) | size(v,2) ~= size(z,2) 
         v = reshape(v,size(z,1),size(z,2));
      end
      z = z + v;
   end
end

d2L = triu(z);
if ~issparse(d2L), d2L = sparse(d2L); end

% MODIFICATION LOG:
%
% 030708 ango Wrote file, based on cpt_d2L.m
% 031204 hkh  Fix to handle vector instead of matrix returned from AD
% 041126 ango Fixed missing addition
% 070823 ango Return sparse triu matrix, not static dense vector representation
% 090813 med  mlint check