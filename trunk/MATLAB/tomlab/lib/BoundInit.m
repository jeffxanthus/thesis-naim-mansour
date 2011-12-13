% function [n, x_k, xkm_1, xEqual, x_L, x_U, Prob] = BoundInit(Prob)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., Sweden. $Release: 7.2.0$
% Written June 28, 1999.    Last modified Jul 20, 2009.

function [n, x_k, x_km1, xEqual, x_L, x_U, Prob] = BoundInit(Prob)

if nargin < 1
   error('BoundInit must have 1 parameter Prob as input');
end

x_k = Prob.x_0(:);  
x_L = Prob.x_L(:);
x_U = Prob.x_U(:);

n = length(x_k);                  % # of variables

if n == 0
   n = max(length(x_L),length(x_U));
end

if isempty(x_k), x_k = zeros(n,1); end

if isempty(x_U),x_U= Inf*ones(n,1); end
if isempty(x_L),x_L=-Inf*ones(n,1); end

if length(x_L) < n
   if length(x_L)==1 
      x_L=x_L*ones(n,1);
   else 
      x_L=[x_L;-Inf*ones(n-length(x_L),1)]; 
   end
end
if length(x_U) < n
   if length(x_U)==1 
      x_U=x_U*ones(n,1);
   else
      x_U=[x_U;Inf*ones(n-length(x_U),1)]; 
   end
end

if any(isnan(x_k)), x_k = zeros(n,1); end

if n==length(x_k)
   x_k=min(x_k,x_U); % Fix that x is below upper bounds
   x_k=max(x_k,x_L); % Fix that x is greater or equal to lower bounds
   % If x in [x_L,x_L+xTol] or [x_U-xTol,x_U], fix x on bounds
end

x_km1    = inf*ones(n,1);
Prob.x_0 = x_k;
Prob.x_L = x_L;
Prob.x_U = x_U;
Prob.N   = n;

xEqual=eq(x_L,x_U);