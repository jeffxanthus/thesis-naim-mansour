function [obj,grad,hess] = linobj( x, Prob )
%        [obj,grad,hess] = linobj( x, Prob )
%        computes the objective value, gradient and diagonal Hessian
%        of a separable convex function, for use with pdco.m.
%        This is an example LINEAR objective function.

%-----------------------------------------------------------------------
% 20 Jun 1997: First test case for pdsco.m.
%              obj  =  c'x  +  1/2 x'Hx   for some c and diagonal H.
%              grad =  c    +  Hx
% 30 Mar 1998: Use linear function to test pdsco.m on LPs.  (H = 0)
% 08 May 1998: c = ones(n,1) gives g - A'y - z = 0 initially.
%              Make c(j) = 0, j even.
% 30 Sep 2002: Modified for pdco.m.
%-----------------------------------------------------------------------
   if nargin < 2
      Prob = [];
   end

   n    = length(x);
   c    = ones(n,1);
   if n>=2, c(2:2:n) = zeros(n/2,1); end
   hess = zeros(n,1);
   obj  = c'*x;
   grad = c;

%-----------------------------------------------------------------------
% End of linobj.m
%-----------------------------------------------------------------------