% FDHess
%
% Numerical approximation of the Hessian matrix
%
% Implementation based on the algorithm FD in Practical Optimization, page 343.
%
% function H = FDHess(x, Prob, gx, varargin)
%
% Sparsity pattern in Prob.HessPattern is used.
% Preprocessed efficient sparsity input in Prob.HessIx is utilized

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Nov 24, 1998.  Last modified July 27, 2011.

function H = FDHess(x, Prob, gx, varargin)

if nargin < 3, gx = []; end

global gTol

x = x(:);
n = Prob.N;

if isempty(gx) 
   gx = nlp_g(x, Prob, varargin{:});
end
Pattern = Prob.HessPattern;
HessIx  = Prob.HessIx;

% HTOL is vector of relative intervals, See Pract. Opt. page 130
HTol    = Prob.GradTolH;
NumDiff = Prob.NumDiff;

if isempty(HTol)
   HTol(1:n) = Prob.optParam.DiffInt;
end
if HTol(1) <= 0
   if ~isempty(gTol)
      % Use the result from the gradient estimation
      HTol      = 10*gTol;
   else
      h = Prob.optParam.DiffInt;
      if h < 0
         HTol = 1E-6*ones(n,1);
      else
         HTol = h*ones(n,1);
      end
   end
else
   HTol=HTol(:);
end
if length(HTol) < n
   HTol = HTol(1)*ones(n,1);
end

% Safeguard against x outside simple bounds
ixx = zeros(n,1);
ixx(1:n)=ixx(1:n)+(x(1:n)+HTol(1:n).*(1+abs(x(1:n))) < Prob.x_L(1:n));
ixx(1:n)=ixx(1:n)+(x(1:n)+HTol(1:n).*(1+abs(x(1:n))) > Prob.x_U(1:n));
% If violating bound, make step in opposite direction, changing sign
ixx = ixx(ixx>0);
HTol(ixx)=-HTol(ixx);

if isempty(Pattern)
   H = zeros(n,n);
else
   [ix,iy,iv]=find(Pattern);
   if isempty(ix)
      % Pattern given, but all zeros. Nothing to calculate.
      H = zeros(n,n);
      return;
   end
   iv = double(iv);
end

if isempty(Pattern) & NumDiff <= 0
   for dim = 1:n
      Prob.FDVar = dim;
      z          = x(dim);
      h          = HTol(dim)*(1+abs(z));
      x(dim)     = z + h;
      gx_ph      = nlp_g(x,Prob, varargin{:});
      x(dim)     = z;
      H(:,dim)   = (gx_ph-gx)/h;
   end   
   H = 0.5*(H+H'); % Make perfectly symmetric
elseif isempty(Pattern)
   % Two levels of differentiation, use f(x) values instead
   f             = nlp_f(x,Prob,varargin{:});
   if abs(f) >= 1E300
      % Singular point, try disturbing x
      ss                   = ones(n,1);
      ss(rand(n,1) <= 0.5) = -1;
      h                    = HTol;
      if isempty(h), h = 1E-8; end
      x                    = x + 0.1*h.*ss;
      f                    = nlp_f(x, Prob, varargin{:});
   end
   % Increase tolerance by factor 10
   HTol = 10*HTol;
   for i = 1:n
      z          = x(i);
      Prob.FDVar = i;
      h          = HTol(i)*(1+abs(z));
      x(i)       = z + h;
      f1         = nlp_f(x,Prob,varargin{:});
      x(i)       = x(i) + h;
      %x(i)       = z - h;
      f11        = nlp_f(x,Prob,varargin{:});
      H(i,i)     = ((f-f1)+(f11-f1))/h^2;
      %H(i,i)     = ((f1-f)-(f-f11))/h^2;
      for j = i+1:n
         x(i)       = z + h;
         v          = x(j);
         Prob.FDVar = [i,j];
         h2         = HTol(j)*(1+abs(v));
         x(j)       = v + h2;
         f12        = nlp_f(x,Prob,varargin{:});
         x(i)       = z;
         Prob.FDVar = j;
         f2         = nlp_f(x,Prob,varargin{:});
         x(j)       = v;
         H(j,i)     = ((f12-f2)+(f-f1))/(h*h2);
         H(i,j)     = H(j,i);
      end
   end
elseif ~isempty(HessIx)
   ixS=[]; iyS=[]; ivS=[];
   mx = max(HessIx);
   for k = 1:mx
      CI          = find(HessIx==k);
      Prob.FDVar  = CI;
      z           = x(CI);
      h           = HTol(CI).*(1+abs(z));
      if NumDiff > 0 % Also the gradient is estimated numerically
         Prob.g_k = gx; 
         for j = 1:length(CI)
             % only estimate NaN elements in nlp_g
             Prob.g_k(ix(find(CI(j)==iy))) = NaN;  
         end
         Prob.cols = find(isnan(Prob.g_k));
      else
         if length(CI) == 1
            Prob.cols = ix(find(CI==iy));
         else
            Prob.cols = ix(find(any([ones(length(iy),1)*CI==iy*ones(1,mx)]')));
         end
      end
      x(CI)       = z + h;
      gx_ph       = nlp_g(x,Prob, varargin{:});
      x(CI)       = z;
      for j = 1:length(CI)
          iz      = find(CI(j)==iy);
          ixz     = ix(iz);
          tmp     = (gx_ph(ixz)-gx(ixz))/h(j);
          iv(iz)  = tmp;
          t       = find(ixz~=CI(j));
          if ~isempty(t)
             ixS  = [ixS;CI(j)*ones(length(t),1)];
             iyS  = [iyS;ixz(t)];
             ivS  = [ivS;tmp(t)];
          end
      end
   end
   H=sparse([ix;ixS],[iy;iyS],[iv;ivS],n,n);
elseif NumDiff > 0
   % Assumption. Only upper rectangle is stored in Pattern
   % Two levels of differentiation, use f(x) values instead
   [ix,iy,iv]=find(Pattern);
   ixS = zeros(length(ix),1); 
   iyS = zeros(length(ix),1); 
   ivS = zeros(length(ix),1);
   % Increase tolerance by factor 10
   HTol = 10*HTol;
   nS = 0;
   f                = nlp_f(x,Prob,varargin{:});
   if n == 1
      ic          = 1;
   else
      ic          = find(any(Pattern));
   end
   for ii = 1:length(ic)
      i          = ic(ii);
      iz         = find(i==iy);
      ixz        = ix(iz);
      m          = length(ixz)-1;
      Prob.FDVar = i;
      z          = x(i);
      h          = HTol(i)*(1+abs(z));
      x(i)       = z + h;
      f1         = nlp_f(x,Prob,varargin{:});
      x(i)       = x(i) + h;
      f11        = nlp_f(x,Prob,varargin{:});
      y          = zeros(m+1,1);
      y(end)     = ((f-f1)+(f11-f1))/h^2;
      for jj = 1:m
         j          = ixz(jj);
         x(i)       = z + h;
         v          = x(j);
         Prob.FDVar = [i,j];
         h2         = HTol(j)*(1+abs(v));
         x(j)       = v + h2;
         f12        = nlp_f(x,Prob,varargin{:});
         x(i)       = z;
         Prob.FDVar = j;
         f2         = nlp_f(x,Prob,varargin{:});
         x(j)       = v;
         y(jj)      = ((f12-f2)+(f-f1))/(h*h2);
      end
      x(i)          = z;
      iv(iz)        = y;
      if m > 0
         ixS(nS+1:nS+m) = i;
         iyS(nS+1:nS+m) = ixz(1:m);
         ivS(nS+1:nS+m) = y(1:m);
         nS = nS + m;
      end
   end   
   H=sparse([ix;ixS(1:nS)],[iy;iyS(1:nS)],[iv;ivS(1:nS)],n,n);
else
   % Assumption. Only upper rectangle is stored in Pattern
   [ix,iy,iv]=find(Pattern);
   ixS=[]; iyS=[]; ivS=[];
   for dim = 1:n
      iz            = find(dim==iy);
      if ~isempty(iz)
         Prob.FDVar = dim;
         z          = x(dim);
         h          = HTol(dim)*(1+abs(z));
         ixz        = ix(iz);
         Prob.cols  = ixz;
         x(dim)     = z + h;
         if NumDiff > 0 % Also the gradient is estimated numerically
            Prob.g_k      = gx; 
            Prob.g_k(ixz) = NaN;  % only estimate NaN elements in nlp_g
         end
         gx_ph      = nlp_g(x,Prob, varargin{:});
         x(dim)     = z;
         tmp        = (gx_ph(ixz)-gx(ixz))/h;
         iv(iz)     = tmp;
         t          = find(ixz~=dim);
         if ~isempty(t)
            ixS=[ixS;dim*ones(length(t),1)];
            iyS=[iyS;ixz(t)];
            ivS=[ivS;tmp(t)];
         end
      end
   end   
   H=sparse([ix;ixS],[iy;iyS],[iv;ivS],n,n);
end

% MODIFICATION LOG
%
% 990306  hkh  Safeguard against x slightly outside bounds
% 000911  hkh  Major revision
% 001022  hkh  Must use n in sparse command to get right dim of H
% 020416  hkh  Use Prob.N for length of x, if x is longer (minimax)
% 030206  ango Check for all zeros in HessPattern
% 040125  hkh  Set double(iv) to avoid treatment at logical array in Matlab 6.5
% 040407  hkh  Send Prob.FDVar telling indicies of perturbed variables
% 040407  hkh  If HessIx defined, use a more efficient method with less calls
% 040407  hkh  Set Prob.g_k with NaN, more efficient two-level differentiation
% 040407  hkh  Assume and utilize only upper triagonal of HessPattern
% 040412  hkh  New algorithm using f(x) if doing differentiation twice
% 040413  hkh  Check length of x_L and x_U
% 041210  hkh  Use 10*gTol, not gTol
% 060218  hkh  Set HTol = 100*HTol for 2 levels of differentiation
% 061117  hkh  More efficient code for safeguard against x outside simple bounds
% 080603  hkh  If f(x) singular, disturb x, recompute f(x) (2 levels of diff)
% 090813  med  mlint check
% 110727  hkh  Use 10*HTol for 2 levels of differentiation
