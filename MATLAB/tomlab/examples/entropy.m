function [obj,grad,hess] = entropy( x, Prob )
%        [obj,grad,hess] = entropy( x, Prob )
%        computes the objective value, gradient and diagonal Hessian
%        of a separable convex function, for use with pdsco.m.
%        This is an example objective function.

%-----------------------------------------------------------------------
% 15 May 1998: Entropy function suggested by Michael Grant.
%              obj     =  sum  x(j)*log x(j)
%              grad(j) =  1 + log x(j)
%              hess(j) =  1 / x(j)
%-----------------------------------------------------------------------

   n    = length(x);
   logx = log(x);

   obj  = sum( x.*logx );
   grad = 1 + logx;
   hess = 1 ./ x;

%-----------------------------------------------------------------------
% End of entropy.m
%-----------------------------------------------------------------------
