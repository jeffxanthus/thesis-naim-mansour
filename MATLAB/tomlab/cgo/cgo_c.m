% function c = cgo_c(X, Prob, varargin)
%
% Compute noncostly constraints c(x) 
% for set of n column vectors x_i in X
%
% OUTPUT
% c     mN x n matrix with costly constraint values 
%       c(j, x_i), j=1,...,mN, i=1,...,n

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written June 30, 2008.   Last modified May 13, 2009.

function c = cgo_c(X, Prob, varargin)

n = size(X,2);
if Prob.simType == 2
   c = zeros(Prob.CGO.mNonLin,n);
   if Prob.CGO.cNARG > 2
      for i=1:n
         c(:,i) = feval(Prob.CGO.c,X(:,i),Prob,varargin{:});
      end
   elseif Prob.CGO.cNARG > 1
      for i=1:n
         c(:,i) = feval(Prob.CGO.c,X(:,i),Prob);
      end
   else
      for i=1:n
         c(:,i) = feval(Prob.CGO.c,X(:,i));
      end
   end
else
   c = zeros(Prob.mNonLin,n);
   if Prob.NARG(4) > 2
      for i=1:n
          c(:,i) = feval(Prob.FUNCS.c,X(:,i),Prob,varargin{:});
      end
   elseif Prob.NARG(4) > 1
      for i=1:n
          c(:,i) = feval(Prob.FUNCS.c,X(:,i),Prob);
      end
   else
      for i=1:n
          c(:,i) = feval(Prob.FUNCS.c,X(:,i));
      end
   end
end

% MODIFICATION LOG
%
% 080630 hkh  Written
% 090513 hkh  Also handle the case userconstraints(x), no Prob argument
