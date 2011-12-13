% function LS = StatLS(x_k, r_k, J_k);
%
% Compute parameter statistics
%
% x_k     Optimal parameter vector, length n
% r_k     Residual vector, length m
% J_k     Jacobian matrix, length m by n
%
% OUTPUT:
% Structure LS with fields
%
% SSQ      Sum of squares: r_k'*r_k;
% Covar    Covariance matrix: Inverse of J' *diag(1./(r_k'*r_k)) * J
% sigma2   Estimate squared standard deviation of problem, 
%          SSQ / Degrees of freedom, i.e. SSQ/(m-n);
% Corr     Correlation matrix: Normalized Covariance matrix
%          Cov./(CovDiag*CovDiag'), where CovDiag = sqrt(diag(Cov));
% StdDev   Estimated standard deviation in parameters: CovDiag*sqrt(sigma2);
% x        =x_k, the input x
% ConfLim  95% Confidence limit (roughly) assuming normal distribution of errors
%          ConfLim = 2*LS.StdDev;
% CoeffVar The coefficients of variation of estimates: StdDev./x_k

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Aug 13, 2002.    Last modified May 11, 2004.

function LS = StatLS(x_k, r_k, J_k)

% Parameter statistics, get Jacobian from Result structure

[m,n] = size(J_k);

mn = m*n;

if mn < 10000
   J=full(J_k); % Convert J from sparse format to full so svd() works
   [u,s,v]=svd(J);
   s = diag(s);
   ix = find(s>0);
   sI = s;
   sI(ix) = 1./s(ix).^2;
   % Covariance matrix: Inverse of J' * J
   lensI = length(sI);
   lenv = size(v,1);
   if lensI < lenv
       sI = [sI; zeros(lenv-lensI,1)];
   end
   Cov = v * diag(sI) * v';
else
   JJ = J_k'*J_k;
   Cov = pinv(JJ);
end

LS.SSQ=r_k'*r_k;

LS.Covar = Cov;

% Estimate squared standard deviation of problem, SSQ / Degrees of freedom
if m > n
   sigma2 = LS.SSQ/(m-n);
else
   sigma2 = LS.SSQ;
end

LS.sigma2 = sigma2;

% Square root of diagonal elements of the covariance matrix

CovDiag = sqrt(diag(Cov));

% Correlation matrix: Normalized Covariance matrix
LS.Corr = Cov./(CovDiag*CovDiag');


% Estimated standard deviation in parameters:
LS.StdDev = CovDiag*sqrt(sigma2);

LS.x = x_k;

% 95% Confidence limit (roughly) assuming normal distribution of errors
LS.ConfLim = 2*LS.StdDev;

% The coefficients of variation of estimates: StdDev./x_k
LS.CoeffVar=LS.StdDev ./ x_k;

% MODIFICATION LOG
%
% 040402  hkh  Add comments about output parameters
% 040511  hkh  Avoid memory overflow for large residual problems