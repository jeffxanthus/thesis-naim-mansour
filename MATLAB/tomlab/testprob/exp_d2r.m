% function d2r = exp_d2r(x, r, J, Prob)
%
% Computes the 2nd part of the second derivative for the
% nonlinear least squares exponential fitting problem at the point x

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.2.0$
% Written Oct 24, 1998.  Last modified Jun 5, 2008.

function d2r = exp_d2r(x, r, J, Prob)

global LS_A

E=LS_A;

x=x(:);
n=length(x);

[SepAlg, p, wType, eType] = expGet1(Prob);

d2r=[];
y  = Prob.LS.y;
t  = Prob.LS.t;

m=size(y,1);

if ~SepAlg
    % Note: b = lambda, a = alpha, c = beta
    alpha=x(p+1:p+p);
    d2r=zeros(n,n);
    if eType==1
        if wType==0 % No weighting
            rW = r .* t;
            ab = -rW' * E;
            bb = (rW.*t)' * ( (ones(m,1) * alpha') .* E);
            d2r = [diag(bb), diag(ab); diag(ab), zeros(p,p)];
        end
    end
end

% MODIFICATION LOG:
%
% 981024  hkh  Fix 2nd derivative for eType==1, wType==0
% 990216  hkh  Expand routine to handle eType=5; distribution estimation