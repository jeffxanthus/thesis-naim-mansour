% function f = ego05_f(x, Prob)
%
% Expected improvement function

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2005 by Tomlab Optimization Inc., $Release: 3.0.0$
% Written Sep 29, 1998.   Last modified Jan 6, 2006.

function f = ego05_f(x, Prob)

invR   = Prob.EGO.invR;
yMin   = Prob.EGO.yMin;
k      = Prob.EGO.k;
p      = Prob.EGO.p;
n      = size(invR,1);

if length(p) > 1
   % p in R^d, same length as x, all p in ]0,2]
   r = exp(-( Prob.EGO.theta'*(abs( ...
       x(1:k)*ones(1,n)-Prob.EGO.X ).^(p*ones(1,n))) ))';
else
   r = exp(-( Prob.EGO.theta'*(abs( ...
       x(1:k)*ones(1,n)-Prob.EGO.X ).^p) ))';
end

%r = zeros(n,1);
%for i = 1:n
%   r(i) = exp(-(  theta'*(abs( x-X(:,i) )).^p ));
%end

Rr = invR * r;

%yHat = my + r'*invR*(Prob.EGO.y-my);

yHat = Prob.EGO.my + Rr'*(Prob.EGO.y-Prob.EGO.my);

tmp = 1-r'*Rr + (1-sum(Rr))^2/sum(sum(invR));

s = sqrt(Prob.EGO.sigma2*max(1E-300,tmp));

%  s=1;  % *********** NOTE ******************

% ML function. Minus sign to make it a minimization problem
if s > 0
   f = -( (yMin - yHat)*phi((yMin-yHat)/s) + s*fi((yMin-yHat)/s) );
else
   f = 1E10;
end
if ~isfinite(f), f=1E10; end

switch(Prob.CGO.EITRANSFORM)
    case 1,
        if(f>0)
            f = 1e10;
        else
            f = -log(-f+1e-5);
        end 
    case 2,
        if(f>0)
            f = 1e10;
        else
            % Is 1e-5 an apropriate marigin to 0?
            f = -1/(f+1e-5);
        end
end

function y = phi(x) % Normal distribution function
y = 0.5  + 0.5 * erf(x/sqrt(2));

function y = fi(x) % Normal density function
y = 1/sqrt(2*pi)*exp(-x^2/2);

% MODIFICATIONS
% 040930  hkh  Safe guard, set f=1E10 if not finite
% 050218  hkh  Clean up
% 060106  frhe Added transformations through EITRANSFORM.