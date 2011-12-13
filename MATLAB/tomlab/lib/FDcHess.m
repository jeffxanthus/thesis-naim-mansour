% FDcHess
%
% Numerical approximation of the nonlinear constraints Hessian
% matrix.
%
% Based on FDHess.m which implements a version of the algorithm FD
% in Practical Optimization, page 343.
%
% function d2c = FDcHess(x,lam,Prob,dc,varargin)
%
% Sparsity pattern in Prob.ConsPattern is used.
% Preprocessed efficient sparsity input in Prob.ConIx is utilized

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Jan 13, 2003.   Last modified July 27, 2011.

function d2c = FDcHess(x,lam,Prob,dc,varargin)

if nargin < 4, dc = []; end

global gTol

x = x(:);
if isempty(dc)
   dc = nlp_dc(x,Prob,varargin{:});
end

lam      = lam(:);
n        = Prob.N;
Pattern  = Prob.ConsPattern;
ConsDiff = Prob.ConsDiff;
HTol     = Prob.GradTolH;

% Step length, check Prob.GradTolH first

% ... then Prob.optParam.DiffInt, or default 1e-6
if isempty(HTol)
   HTol(1:n) = Prob.optParam.DiffInt;
end
if HTol(1) <= 0
   if ~isempty(gTol)
      % Use the result from the gradient estimation
      HTol      = gTol;
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
% Fix length, using first element if not n elements
if length(HTol) < n
   HTol = HTol(1)*ones(n,1);
end

% Safeguard against x outside simple bounds
% (DANGEROUS for central differences, check for 2*step)
ixx = zeros(n,1);
ixx(1:n)=ixx(1:n)+(x(1:n)+2*HTol(1:n).*(1+abs(x(1:n))) < Prob.x_L(1:n));
ixx(1:n)=ixx(1:n)+(x(1:n)+2*HTol(1:n).*(1+abs(x(1:n))) > Prob.x_U(1:n));
% If violating bound, make step in opposite direction, changing sign
ixx = ixx(ixx>0);
HTol(ixx)=-HTol(ixx);

d2c = sparse(n,n);

if isempty(Pattern) & ConsDiff <= 0
   % Case 1: One level of differentiation, no pattern
   ix = find(lam'~=0);
   if ~isempty(ix)
      Prob.rows = ix;
      for i = 1:n
        z          = x(i);
        Prob.FDVar = i;
        h          = HTol(i)*(1+abs(z));
        x(i)     = z + h;
        dc_ph      = nlp_dc(x,Prob,varargin{:});
        x(i)     = z;
        % Consider only constraints for which lam is nonzero
        for k=ix
          d2c(:,i) = d2c(:,i)+ (lam(k)/h)*( dc_ph(k,:)-dc(k,:) )';
        end
      end
   end
   % Make symmetric matrix
   d2c=0.5*(d2c+d2c');
elseif isempty(Pattern)
   % Case 2: Two levels of differentiation, use c(x) values instead, no pattern
   %HTol = 100*HTol;
   HTol  = 5*HTol;
   ix    = find(lam'~=0);
   if ~isempty(ix)
      Prob.rows    = ix;
      c            = nlp_c(x,Prob,varargin{:});
      for i = 1:n
        z          = x(i);
        Prob.FDVar = i;
        h          = HTol(i)*(1+abs(z));
        x(i)       = z + h;
        c1         = nlp_c(x,Prob,varargin{:});
        x(i)       = x(i) + h;
        c11        = nlp_c(x,Prob,varargin{:});
        d2c(i,i) = sum((lam(ix)/h^2).*((c(ix)-c1(ix))+(c11(ix)-c1(ix))));
        %for k=ix
        %  d2c(i,i) = d2c(i,i)+(lam(k)/h^2)*((c(k)-c1(k))+(c11(k)-c1(k)));
        %end
        for j = i+1:n
           x(i)       = z + h;
           v          = x(j);
           Prob.FDVar = [i,j];
           h2         = HTol(j)*(1+abs(v));
           x(j)       = v + h2;
           c12        = nlp_c(x,Prob,varargin{:});
           x(i)       = z;
           Prob.FDVar = j;
           c2         = nlp_c(x,Prob,varargin{:});
           x(j)       = v;
           hh         = h*h2;
           % Consider only constraints for which lam is nonzero
           d2c(j,i)   = sum((lam(ix)/hh).*((c12(ix)-c2(ix))+(c(ix)-c1(ix))));
           %for k=ix
           %  d2c(j,i) = d2c(j,i)+(lam(k)/hh)*((c12(k)-c2(k))+(c(k)-c1(k)));
           %end
           d2c(i,j)   = d2c(j,i);
        end
      end
   end
elseif ~isempty(Pattern) & ConsDiff > 0
   % Case 3: Two levels of differentiation, use c(x) values instead
   % Nonempty Prob.ConsPattern... 
   %HTol = 100*HTol;
   HTol  = 5*HTol;
   ix    = find(lam'~=0);
   if ~isempty(ix)
      if length(ix) == 1
         ic          = find(Pattern(ix,:));
      else
         ic          = find(any(Pattern(ix,:)));
      end
      Prob.rows      = ix;
      c              = nlp_c(x,Prob,varargin{:});
      Prob.cols      = ic;
      for ii = 1:length(ic)
         i           = ic(ii);
         z           = x(i);
         Prob.FDVar  = i;
         h           = HTol(i)*(1+abs(z));
         x(i)        = z + h;
         c1          = nlp_c(x,Prob,varargin{:});
         x(i)        = x(i) + h;
         c11         = nlp_c(x,Prob,varargin{:});
         d2c(i,i)    = sum((lam(ix)/h^2).*((c(ix)-c1(ix))+(c11(ix)-c1(ix))));
         %for k = ix
         %   d2c(i,i) = d2c(i,i)+(lam(k)/h^2)*((c(k)-c1(k))+(c11(k)-c1(k)));
         %end
         for jj = ii+1:length(ic)
            j          = ic(jj);
            x(i)       = z + h;
            v          = x(j);
            Prob.FDVar = [i,j];
            h2         = HTol(j)*(1+abs(v));
            x(j)       = v + h2;
            c12        = nlp_c(x,Prob,varargin{:});
            x(i)       = z;
            Prob.FDVar = j;
            c2         = nlp_c(x,Prob,varargin{:});
            x(j)       = v;
            hh         = h*h2;
            % Consider only constraints for which lam is nonzero
            d2c(j,i)   = sum((lam(ix)/hh).*((c12(ix)-c2(ix))+(c(ix)-c1(ix))));
            %for k = ix
            %   d2c(j,i) = d2c(j,i)+(lam(k)/hh)*((c12(k)-c2(k))+(c(k)-c1(k)));
            %end
            d2c(i,j)   = d2c(j,i);
         end
      end
   end
else
   % Case 4: Nonempty Prob.ConsPattern, one level of differentiation
   ix = find(lam'~=0);
   if ~isempty(ix)
      if length(ix) == 1
         ic = find(Pattern(ix,:));
      else
         ic = find(any(Pattern(ix,:)));
      end
      Prob.rows     = ix;
      Prob.cols     = ic;
      for i = ic
         z          = x(i);
         Prob.FDVar = i;
         h          = HTol(i)*(1+abs(z));
         x(i)       = z + h;
         dc_ph      = nlp_dc(x,Prob,varargin{:});
         x(i)       = z;
         % Consider only constraints for which lam is nonzero
         for k = ix
            d2c(ic,i) = d2c(ic,i)+ (lam(k)/h)*( dc_ph(k,ic)-dc(k,ic) )';
         end
      end
   end
   % Make symmetric matrix
   d2c=0.5*(d2c+d2c');
end

% MODIFICATION LOG
%
% 030113 ango Wrote file
% 030113 hkh  Use Prob.ConsPattern to find which variables to use
% 030113 hkh  Define d2c sparse
% 030114 ango Fix error in x+HTol < Prob.x_L (was -, not +)
% 030128 hkh  Compute find(Lam') outside loop
% 030220 hkh  If ix empty, avoid all computations
% 040408 hkh  Only consider rows in Pattern for nonzero lam
% 040412 hkh  New algorithm using f(x) if doing differentiation twice
% 061117 hkh  More efficient code for safeguard against x outside simple bounds
% 090813 med  mlint check
% 110727 hkh  Incorrect test for Case 3, never executed
% 110727 hkh  Use 5*HTol for 2 levels of differentiation
% 110727 hkh  Avoid loop over ix for single d2c elements, when 2 level of diff

