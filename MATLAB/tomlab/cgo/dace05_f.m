% function f = dace05_f(x, Prob)
%
% Likelihood function

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Sept 27, 1998. Last modified August 24, 2009.

function f = dace05_f(x, Prob)

% Compute the correlation matrix R

d = Prob.dDim;

if length(x) > d
   R=tomsol(28,x(1:d),Prob.DACE.X,x(d+1:2*d));
else
   R=tomsol(28,x,Prob.DACE.X,2);
end

% Old code to get R
% R = getR(x,Prob)

if R(1,1) > 1E99
   % No inverse possible if R(1,1) = 1E100
   f  = 1E10;
   return;
end

n = size(R,1);

% [invR, detR, pRank] = getinvR(R, Prob.optParam.eps_Rank);
[invR, detR, pRank] = tomsol(12, R, Prob.optParam.eps_Rank);

if detR <= 1E-300 | pRank < n
   f = 1E10;
   return;
end

ssinvR = sum(sum(invR));

if ssinvR > 0
   y  = Prob.DACE.y(:);
   my = sum(invR*y)/ssinvR;
   s2 = ((y-my)'*invR*(y-my))/n;
else
   %my = sum(invR*y)/ssinvR;
   %s2 = ((y-my)'*invR*(y-my))/n;
   %if isinf(s2)
   %   keyboard
   %end
        
   f  = 1E10;
   return;
end

% Formula is inverted to create a minimization problem

%f = log( (s2)^(n/2)*sqrt(detR) );

%f = - log( exp(-n/2) / ( (2*pi*s2)^(n/2)*sqrt(det(R)) ) );

f = n/2*log(max(1E-300,s2)) + 0.5*log(max(1E-300,detR));

% --------------------------------------------------------------
function [invR, detR, pRank] = getinvR(R, epsRank)
% --------------------------------------------------------------
% NOT USED NOW
n          = size(R,1);
[U S V]    = svd(R);
S_inv      = zeros(n,1);
S11        = S(1,1);
detR       = S11;
S_inv(1) = 1/S11;
for i = 2:n
    Sii = S(i,i);
    if Sii > epsRank*S11
       pRank    = i;
       detR     = detR*Sii;
       S_inv(i) = 1/Sii;
    else
       break;
    end
end
invR = V(:,1:pRank) * diag(S_inv(1:pRank)) * U(:,1:pRank)';
% sum(invR-inv(R))
% fprintf('invR: pRank %d n %d detR %20.4e\n',pRank,n,detR);

% --------------------------------------------------------------
function R = getR(x,Prob)
% --------------------------------------------------------------
% Old code to get R

X      = Prob.DACE.X;
[k,n]  = size(X);
R      = eye(n);

% If optimizing in log(x), then transform first
expx = exp(x(1:k));

if length(x) > k
   p = x(k+1:2*k);
   for i = 1:n-1
      % x transformed, so log(x) are input
      R(i,i+1:n) = exp(-(  expx'*(...
         abs( X(:,i)*ones(1,n-i)-X(:,i+1:n) ).^(p*ones(1,n-i)))  ));
      R(i+1:n,i) = R(i,i+1:n)';

   end
else
   % p = 2; % NOT USED 
   for i = 1:n-1
      % x transformed, so log(x) are input
      R(i,i+1:n) = exp(-(  expx'*(abs( X(:,i)*ones(1,n-i)-X(:,i+1:n) ).^2)  ));
      R(i+1:n,i) = R(i,i+1:n)';

      %for j = i+1:n
      %   % x not transformed. Then x >= 0 must be used
      %   % R(i,j) = exp(-( x(1:k)'*( abs(X(:,i)-X(:,j)).^p)  ));
      %   
      %   % If x transformed, so log(x) are input
      %   R(i,j) = exp(-(  expx'*(abs( X(:,i)-X(:,j) ).^p)  ));
      %   R(j,i) = R(i,j);
      %end
   end
end

% MODIFICATION LOG
%
% 040309 hkh Revise, tomsol(28) returns R(1,1) 1E100 if no inverse possible
% 050218 hkh Revise, Use SVD to get inv(R) in safe way
% 090824 hkh Minor mlint revision
