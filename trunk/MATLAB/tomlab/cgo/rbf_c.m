% function f = rbf_c(x, Prob)
%
% If Prob.SCALE = 1, then
% rbf_c does rescaling of x before calling the constraint routine
%
% It then computes the distance constraints with the original x
%    ||x-x_i||_2^2 >= (beta*maxD)^2, for selected i in Prob.ixDist


% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 3.0.0$
% Written Mar 28, 2002.   Last modified May 9, 2005.

function c = rbf_c(x, Prob)

n  = Prob.dDim;
N  = Prob.N;
m  = Prob.mNonLin;
ix = Prob.ixDist;
m2 = length(ix);
if isfield(Prob,'MY')
   m3 = 1;
else
   m3 = 0;
end
m1 = m-m2-m3;
c  = zeros(m,1);

if m1 > 0
   if Prob.SCALE
      O = tomsol(9, Prob.xL(1:n), x(1:n) ,Prob.xD(1:n)); 
      if Prob.cNargin>=2
         c(1:m1)=feval(Prob.c, O, Prob);
      else
         c(1:m1)=feval(Prob.c, O);
      end
   else
      if Prob.cNargin>=2
         c(1:m1)=feval(Prob.c, x(1:n), Prob);
      else
         c(1:m1)=feval(Prob.c, x(1:n));
      end
   end
end

% Also handle the following type of problem
% min -x(n+1) 
% s/t ||x-x_i||_2^2-x(n+1)^2 >= 0, 
%     x_L(1:n) <= x(1:n) <= x_U(1:n), x(n+1) >= 0
if m2 > 0
   if N > n
      c(m1+1:m-m3) = (sum((x(1:n)*ones(1,m2)-Prob.user.X(:,ix)).^2))'-x(end)^2;
   else
      c(m1+1:m-m3) = (sum((x*ones(1,m2)-Prob.user.X(:,ix)).^2))';
   end
end
if m3 > 0
   % c(end) = tomsol(26,x(1:n),Prob.CGO.fnStar);
   c(end) = tomsol(26,x(1:n),0);
end

% MODIFICATION LOG:
%
% 020328  hkh  Written
% 030514  hkh  Only use original length of x when calling constraints
% 030514  hkh  Also use this routine when no SCALING
% 040104  hkh  dDim holds original size of problem, n not OK set 
% 050424  hkh  Added distance to sample points constraints
% 050506  hkh  Make to work for expanded problem
% 050509  hkh  Add constraint for my (new RBF coefficient)

