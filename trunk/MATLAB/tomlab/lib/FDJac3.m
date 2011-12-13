% FDJac3
%
% Compute the Jacobian of the residuals in a nonlinear least squares problem
% or the constraint Jacobian for the constraint vector c(x).
%
% function [J,rx] = FDJac3(x, Prob, rxFunc, rx, varargin)
%
% rx = r(x) or c(x), i.e. the current residual or constraint vector,
% if available (otherwise empty)
%
% Sparsity pattern in Prob.JacPattern or Prob.ConsPattern is used.
% Preprocessed efficient sparsity input in Prob.ConIx is utilized
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

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Jan 19, 1999.  Last modified Aug 13, 2009.

function [J,rx] = FDJac3(x, Prob, rxFunc, rx, varargin)

if nargin < 4
   rx = [];
   if nargin < 3
      rxFunc = 'nlp_r';
   end
end

x = x(:);

if strcmp(rxFunc,'nlp_r')
   NLP=1;
elseif strcmp(rxFunc,'nlp_c')
   NLP=2;
else
   NLP=0;
end

if NLP==1 | NLP==0
   if isempty(Prob.JacPattern)
      Pattern=[];
      ConIx = [];
   else
      Pattern=Prob.JacPattern;
      ConIx = Prob.JacIx;
   end
elseif NLP==2
   if isempty(Prob.ConsPattern)
      Pattern=[];
      ConIx = [];
   else
      Pattern=Prob.ConsPattern;
      ConIx = Prob.ConIx;
   end
end

if isempty(rx)
   if NLP==1
      rx = nlp_r(x , Prob, varargin{:} );
   elseif NLP==2
      rx = nlp_c(x, Prob, varargin{:});
   else
      rx = feval(rxFunc, x, Prob, varargin{:});
   end
   ConIx = [];
end

if isempty(rx)
   % If there are no nonlinear constraints then dc should be [].
   J = []; return;
end
   
m = length(rx);
n = Prob.N;

if isempty(Pattern)
   if Prob.LargeScale
      J = sparse(m,n);
   else
      J = zeros(m,n);
   end
else
   [ix,iy,iv]=find(Pattern);
   iv = double(iv);
end

h = 1e-20;

if NLP==1
   if isempty(Pattern)
      for dim = 1:n
         y          = x(dim);
         Prob.FDVar = dim;
         Prob.cols  = dim;
         hS         = h*max(1,abs(y));
         x(dim)     = y + i*hS;
         J(:,dim)   = imag( nlp_r(x,Prob, varargin{:}) )/hS;
         x(dim)     = y;
      end
   elseif 0 && ~isempty(ConIx)
      mx = max(ConIx);
      for k = 1:mx
         CI         = find(ConIx==k);
         nCI        = length(CI);
         Prob.FDVar = CI;
         Prob.cols  = CI;
         if length(CI) == 1
            Prob.rows = ix(find(CI==iy));
         else
            Prob.rows = ix(find(any([ones(length(iy),1)*CI==iy*ones(1,nCI)]')));
         end
         y          = x(CI);
         hS         = h*max(1,abs(y));
         x(CI)      = y + i*hS;
         z          = imag( nlp_r(x,Prob, varargin{:}) )/hS;
         x(CI )     = y;
         for j = 1:length(CI)
             iz = find(CI(j)==iy);
             iv(iz) = z(ix(iz));
         end
      end
   else
      for dim = 1:n
          iz = find(dim==iy);
          if ~isempty(iz)
             ixz        = ix(iz);
             Prob.rows  = ixz;
             y          = x(dim);
             Prob.FDVar = dim;
             Prob.cols  = dim;
             hS         = h*max(1,abs(y));
             x(dim)     = y + i*hS;
             z          = imag( nlp_r(x,Prob, varargin{:}) )/hS;
             x(dim)     = y;
             iv(iz)     = z(ixz);
          end
      end
   end
elseif NLP==2
   if isempty(Pattern)
      for dim = 1:n
         y          = x(dim);
         Prob.FDVar = dim;
         Prob.cols  = dim;
         hS         = h*max(1,abs(y));
         x(dim)     = y + i*hS;
         J(:,dim)   = imag( nlp_c(x,Prob, varargin{:}) )/hS;
         x(dim)     = y;
      end
   elseif 0 && ~isempty(ConIx)
      mx = full(max(ConIx));
      for k = 1:mx
         CI         = find(ConIx==k);
         Prob.FDVar = CI;
         Prob.cols  = CI;
         if length(CI) == 1
            Prob.rows = ix(find(CI==iy));
         else
            Prob.rows = ix(find(any([ones(length(iy),1)*CI==iy*ones(1,mx)]')));
         end
         y          = x(CI);
         hS         = h*max(1,abs(y));
         x(CI)      = y + i*hS;
         z          = imag( nlp_c(x,Prob, varargin{:}) )/hS;
         x(CI )     = y;
         for j = 1:length(CI)
             iz = find(CI(j)==iy);
             iv(iz) = z(ix(iz));
         end
      end
   else
      for dim = 1:n
          iz = find(dim==iy);
          if ~isempty(iz)
             ixz        = ix(iz);
             Prob.rows  = ixz;
             y          = x(dim);
             Prob.FDVar = dim;
             Prob.cols  = dim;
             hS         = h*max(1,abs(y));
             x(dim)     = y + i*hS;
             z          = imag( nlp_c(x,Prob, varargin{:}) )/hS;
             x(dim)     = y;
             iv(iz)     = z(ixz);
          end
      end
   end
else
   if isempty(Pattern)
      for dim = 1:n
          y          = x(dim);
          Prob.FDVar = dim;
          Prob.cols  = dim;
          hS         = h*max(1,abs(y));
          x(dim)     = y + i*hS;
          J(:,dim)   = imag( feval(rxFunc,x,Prob, varargin{:}) )/hS;
          x(dim)     = y;
      end
   elseif 0 && ~isempty(ConIx)
      mx = full(max(ConIx));
      for k = 1:mx
         CI         = find(ConIx==k);
         y          = x(CI);
         Prob.FDVar = CI;
         Prob.cols  = CI;
         if length(CI) == 1
            Prob.rows = ix(find(CI==iy));
         else
            Prob.rows = ix(find(any([ones(length(iy),1)*CI==iy*ones(1,mx)]')));
         end
         hS         = h*max(1,abs(y));
         x(CI)      = y + i*hS;
         z          = imag( feval(rxFunc,x,Prob, varargin{:}) )/hS;
         x(CI )     = y;
         for j = 1:length(CI)
             iz = find(CI(j)==iy);
             iv(iz) = z(ix(iz));
         end
      end
   else
      for dim = 1:n
          iz        = find(dim==iy);
          if ~isempty(iz)
             ixz        = ix(iz);
             Prob.rows  = ixz;
             y          = x(dim);
             Prob.FDVar = dim;
             Prob.cols  = dim;
             hS         = h*max(1,abs(y));
             x(dim)     = y + i*hS;
             z          = imag( feval(rxFunc,x,Prob, varargin{:}) )/hS;
             x(dim)     = y;
             iv(iz)     = z(ixz);
          end
      end
   end
end

if ~isempty(Pattern)
   J=sparse(ix,iy,real(iv),m,n);
else
   J=real(J);
end

% MODIFICATION LOG
%
% 990306  hkh  Safeguard against x slightly outside bounds
% 000910  hkh  Written efficiently, avoid arrays and zero mult. Use 1E-20
% 020416  hkh  Use Prob.N for length of x, if x is longer (minimax)
% 031214  hkh  Use J = sparse(m.n) if Prob.LargeScale = 1
% 040125  hkh  Set double(iv) to avoid treatment at logical array in Matlab 6.5
% 040331  hkh  Return rx as second parameter
% 040406  hkh  If ConIx defined, use a more efficient method with less calls
% 040406  hkh  Send global index FDVar telling index of perturbed variable
% 040406  hkh  Send Prob.FDVar telling indicies of perturbed variables
% 040410  hkh  Send Prob.cols and Prob.rows
% 040505  hkh  Incorrect use of mx as length of CI
% 050124  hkh  ConIx=Prob.JacIx must be set for NLLS problems