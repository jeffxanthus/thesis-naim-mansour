% function X = daceInts(X,x_D,x_U,IV)
%
% daceInts generates integer points in X for the variable indices defined
% in IV. The integer points are generated with equal probability.
%
%
% Input parameters:
% X             Matrix of already sampled points.
% x_D           x_D = x_U-x_L, where x_U,x_L is the upper and lower bounds on x
% x_U           Upper bounds for the variables x.
% IV            Indices for integer variales among the 1,...d variables
%
% Output parameters:
% X             Matrix of sampled points, X(IV,:) integer values between bounds

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Feb 2, 2005.    Last modified Sep 18, 2009.

% ====================================================================
function X = daceInts(X,x_D,x_U,IV)
% ====================================================================
% With equal probability, select integer point, remove duplicates
% IV = IntVars
[d,n]   = size(X);
ix      = ones(n,1);
XM      = (x_D(IV)+1)./max(0,x_D(IV));
X(IV,1) = min(x_U(IV),floor(XM.*X(IV,1)));

for i = 2:n
    X(IV,i) = min(x_U(IV),floor(XM.*X(IV,i)));
    if any(sum(X(:,1:i-1)==X(:,i)*ones(1,i-1))==d)
      ix(i) = 0; % Mark duplicate point
    end
end
ix = find(ix);
if length(ix) < n
   % Only use unique points
   X = X(:,ix);
end

% MODIFICATION LOG:
% 050218  hkh  New function daceInts, round experimental design points
% 050218  hkh  Add check on duplicate points
% 050218  hkh  Use daceInts instead of sampleInts for DACE strategies
% 050322  hkh  Wrong check for duplicates in daceInts, also round only IntVars
% 090918  hkh  Generate with equal probability for every integer value
% 090918  hkh  Make daceInts separate file, avoid duplicates
