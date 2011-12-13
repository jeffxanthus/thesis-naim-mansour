% TOMLAB gateway routine
% Callback from CONOPT
%
% nlp_d2L computes the Hessian
%
%       d2f(x) + lam * d2c(x)
%
% to the Lagrangian function,
%
%   L(x,lam) =   f(x) + lam' * c(x)
%
% function d2L=cpt_d2L(x, lam, Prob)

% Anders Göran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Nov 11, 2002.   Last modified Jan 27, 2003.

function d2L=cpt_d2L(x, lam, Prob)

x=x(:);

z = nlp_H(x,Prob);

if ~isempty(lam)
    m2    = Prob.m2;
    l     = lam(1:m2);             % Original constraints multipliers
    ix    = Prob.cexidx(m2+1:end); % The original indices of the extra constraints
    l(ix) = max(abs(lam(ix)),abs(lam(m2+1:end)));
    if isempty(z)
        z = nlp_d2c(x,l,Prob);
    else
        z = z + nlp_d2c(x,l,Prob);
    end
end

% A dense vector with static sparse information used
% z = sparse(tril(z));
z   = tril(z);
z   = z(:);
d2L = full( z(Prob.D2Lidx) );

% MODIFICATION LOG:
%
% 021111 ango Wrote file
% 021112 ango Added handling of sparse Hessians
% 021223 hkh  Faster to avoid tests before doing sparse/full/transpose
% 021229 hkh  Change comments
% 021231 ango Always return sparse Hessian and transpose removed
% 030127 hkh  Display warning if not symmetric
% 030604 ango Re-wrote, based on nlp_d2L