% function ix = findpatt(C)
%
% findpatt computes the sequence of calls needed to approximate the Jacobian
% of the constraints, where the dependencies are marked as 1s in the input
% matrix C
%
% INPUT:
% C    A 0-1 m by n-matrix, sparse or dense, with the constraint pattern
%      see the description of Prob.ConsPattern
%
% OUTPUT:
% ix   an n-vector, with indices 1,2,...,mx, where mx are the number of calls
%      needed to compute the full numerical Jacobian.
%      -Inf indicates a zero column in C for this variables
%      An ix-value i means that the variables are independent with respect to
%      the constraints they influence, and therefore the gradients may be
%      computed simultaneously

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Apr 5, 2004.  Last modified Apr 9, 2004.

function ix = findpatt(C)

j     = 0;
[m,n] = size(C);
if n > 200
   C = sparse(C);
end
if m == 1
   ix = full(double(C==0));
else
   ix = full(double(all(C==0)));
end

ix(find(ix)) = -Inf;

for i=1:n
   if ix(i) == 0
      j = j+1;
      ix(i) = j;
      r = C(:,i);
      y = find(r);
      while ~isempty(y)
         if length(y) == 1
            z = C(y,i+1:end);
         else
            z = any(C(y,i+1:end));
         end
         v = find(z==0 & ix(i+1:end)==0);
         %v = find(z+ix(i+1:end)==0); SLOWER
         if isempty(v), break; end
         k     = i+v(1);
         ix(k) = j;
         if all(ix), break; end
         r = r | C(:,k);
         y = find(r);
      end
   end
end

% MODIFICATION LOG
%
% 040405  hkh  Algorithm formulated and written
% 040406  hkh  Clean-up and comments
% 040408  hkh  Special case if only one row
% 040409  hkh  Make ix full