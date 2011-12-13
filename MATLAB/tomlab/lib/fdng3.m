% Finite-Difference Numerical approximation of Gradient. (FDNG)
%
% function g = fdng3(x, Prob, g, fx, varargin)
%
% Implementation based on the paper:
% "USING COMPLEX VARIABLES TO ESTIMATE DERIVATIVES OF REAL FUNCTIONS",
% William Squire, George Trapp, SIAM Review, Vol. 10, No. 1, 1998, pp. 100-112
%
% For the coding of the TOMLAB m-files, the following paper has good hints:
%
% Martins, Kroo, Alonso:
% An Automated Method for Sensitivity Analysis using Complex Variables
% 38th Aerospace Sciences Meeting and Exhibit, January 10-13, 2000, Reno, NV
%
% See the papers above before coding your own functions to use this type of
% derivative strategy
%
% If g is nonempty, estimate any elements of g set to NaN.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Jan 19, 1999.  Last modified Apr 11, 2004.

function g = fdng3(x, Prob, g, fx, varargin)

n = Prob.N;

if isempty(g)
   ALLg = 1;
   g = zeros(n,1);
else
   ALLg = 0;
   ix = find(isnan(g));
   if isempty(ix), return, end
end

if nargin < 4, fx = []; end

x     = x(:);

if isempty(fx)
   fx = nlp_f(x,Prob, varargin{:});
end

h = 1e-20;

if ALLg
   for dim = 1:n
      Prob.FDVar = dim;
      z          = x(dim);
      hS         = h*max(1,abs(z));
      x(dim)     = z + i*hS;
      g(dim)     = imag( nlp_f(x,Prob, varargin{:}) )/hS;
      x(dim)     = z;
   end
else
   for k = 1:length(ix)
      dim        = ix(k);
      Prob.FDVar = dim;
      z          = x(dim);
      hS         = h*max(1,abs(z));
      x(dim)     = z + i*hS;
      g(dim)     = imag( nlp_f(x,Prob, varargin{:}) )/hS;
      x(dim)     = z;
   end
end

g=real(g);

% MODIFICATION LOG
%
% 990306  hkh  Safeguard against x slightly outside bounds
% 990626  hkh  Avoid feval
% 000911  hkh  Clean up, written more efficiently
% 020416  hkh  Use Prob.N for length of x, if x is longer (minimax)
% 030127  hkh  Added g as input, estimate NaN elements
% 030131  hkh  Cannot use i for loop index, change to k
% 040411  hkh  Send Prob.FDVar with the variable index changed