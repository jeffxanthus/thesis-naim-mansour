% function [f,c] = cgo_fc(X, Prob, varargin)
%
% Compute costly function f(x) 
% (and costly constraint mC-vector Cc(x) if defined)
% for set of n column vectors x_i in X
% Using parfor loops, utilized in Parallel Computing Toolbox
%
% OUTPUT
% f     n costly function values f(x_i), i=1,...,n
% c     mC x n matrix with costly constraint values 
%       c(j, x_i), j=1,...,mC, i=1,...,n

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written June 29, 2008.   Last modified July 23, 2011.

function [f,c] = cgo_fc(X, Prob, varargin)

n   = size(X,2);
%if n > 1 & all(Prob.MATLAB >= [7,5])
if n > 1 & Prob.Threads > 1
   [f,c] = cgo_fc_par(X, Prob, varargin);
else
   f = zeros(n,1);

   if Prob.simType > 0
      c    = zeros(Prob.mNonLin,n);
      if Prob.NARG(10) > 2
         for i=1:n
            [f(i),c(:,i)] = feval(Prob.FUNCS.fc,X(:,i),Prob,varargin{:});
         end
      elseif Prob.NARG(10) > 1
         for i=1:n
            [f(i),c(:,i)] = feval(Prob.FUNCS.fc,X(:,i),Prob);
         end
      else
         for i=1:n
            [f(i),c(:,i)] = feval(Prob.FUNCS.fc,X(:,i));
         end
      end
   else
      c     = [];
      if Prob.NARG(1) > 2
         for i=1:n
            f(i)       = feval(Prob.FUNCS.f,X(:,i),Prob,varargin{:});
         end
      elseif Prob.NARG(1) > 1
         for i=1:n
            f(i)       = feval(Prob.FUNCS.f,X(:,i),Prob);
         end
      else
         for i=1:n
            f(i)       = feval(Prob.FUNCS.f,X(:,i));
         end
      end
   end
end
   


% MODIFICATION LOG
%
% 080629 hkh  Written
% 090513 hkh  Also handle case userfunc(x) with no Prob
% 090923 hkh  Change for to parfor loops
% 091001 hkh  Use parfor only for new Matlab versions, test on Prob.MATLAB
% 091021 ango Break out parfor code to avoid problems with older Matlab versions
% 100910 ango Matlab version check changed
% 110723 hkh  Only use parfor if loop items > 1
