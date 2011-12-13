% function f = rbf_d2c(x, Prob)
%
% rbf_d2c does rescaling of x before calling the second order constraint
% routine, if Prob.SCALE = 1
%
% It then computes the constraint Jacobian of the distance constraints 
% with the original x
%    ||x-x_i||_2^2 >= (beta*maxD)^2, for selected i in Prob.ixDist

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 3.0.0$
% Written Mar 28, 2002.   Last modified May 9, 2005.

function d2c = rbf_d2c(x, lam, Prob)

N   = Prob.N;
n   = Prob.dDim;
m   = Prob.mNonLin;
m2  = length(Prob.ixDist);
m1  = m-m2;
d2c = zeros(N,N);

if m1 > 0
   if Prob.SCALE
      O  = tomsol(9, Prob.xL(1:n), x(1:n) ,Prob.xD(1:n)); 
      if Prob.d2cNargin>=3
         d2c(1:n,1:n)=feval(Prob.d2c, O, lam(1:m1), Prob);
      elseif Prob.d2cNargin==2
         d2c(1:n,1:n)=feval(Prob.d2c, O, lam(1:m1)); 
      end
   else
      if Prob.d2cNargin>=3
         d2c(1:n,1:n)=feval(Prob.d2c, x(1:n), lam(1:m1), Prob);
      elseif Prob.d2cNargin==2
         d2c(1:n,1:n)=feval(Prob.d2c, x(1:n), lam(1:m1)); 
      end
   end
end

% Also handle the following type of problem
% min -x(n+1) 
% s/t ||x-x_i||_2^2-x(n+1)^2 >= 0, 
%     x_L(1:n) <= x(1:n) <= x_U(1:n), x(n+1) >= 0

if m2 > 0
   if N > n
      V      = eye(N,N);
      V(N,N) = -1;
      d2c    = d2c + 2*sum(lam(m1+1:m1+m2))*V;
   else
      d2c    = d2c + 2*sum(lam(m1+1:m1+m2))*eye(N,N);
   end
end
if m3 > 0
   % Skip this right now
   % c(end) = tomsol(26,y,Prob.CGO.fnStar);
end

   
% MODIFICATION LOG:
%
% 030514  hkh  Only use original length of x when calling constraints
% 030514  hkh  Also use this routine when no SCALING
% 050424  hkh  Added d2c for distance to sample point constraints
% 050506  hkh  Make to work for expanded problem
% 050509  hkh  Add constraint for my (new RBF coefficient)
