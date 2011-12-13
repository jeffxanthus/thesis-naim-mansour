% FDrHess
%
% Numerical approximation of the residuals Hessian matrix.
%
%   sum(i=1:m) r_i*d2r_i
%
% Based on FDcHess.m which is based on FDHess.m which implements
% a version of the algorithm FD in Practical Optimization, page 343.
%
% function d2r = FDrHess(x,Prob,r,J,varargin)
%
% Sparsity pattern in Prob.JacPattern is used.
% NOT preprocessed efficient sparsity input in Prob.ConIx is utilized

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Jan 21, 2005.   Last modified Aug 13, 2009.

function d2r = FDrHess(x,Prob,r,J,varargin)

global gTol

if nargin < 5
  J = [];
  if nargin < 4
    r = [];
  end
end

x = x(:);

if isempty(J)
  J = nlp_J(x,Prob,varargin{:});
end
if isempty(r)
  r = nlp_r(x,Prob,varargin{:});
end
  
r        = r(:);
n        = Prob.N;
Pattern  = Prob.JacPattern;
NumDiff  = Prob.NumDiff;
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

d2r = sparse(n,n);

if isempty(Pattern) & NumDiff <= 0
   ix = find(r');
   if ~isempty(ix)
      Prob.rows = ix;
      for i = 1:n
        z          = x(i);
        Prob.FDVar = i;
        h          = HTol(i)*(1+abs(z));
        x(i)       = z + h;
        J_ph       = nlp_J(x,Prob,varargin{:});
        x(i)       = z;
    
        % Consider only residuals for which r is nonzero
        for k=ix
          d2r(:,i) = d2r(:,i)+ (r(k)/h)*( J_ph(k,:)-J(k,:) )';
        end
      end
   end
   % Make symmetric matrix
   d2r=0.5*(d2r+d2r');
elseif isempty(Pattern)
   % Two levels of differentiation, use r(x) values instead
   ix = find(r');
   if ~isempty(ix)
      Prob.rows    = ix;
      for i = 1:n
        z          = x(i);
        Prob.FDVar = i;
        h          = HTol(i)*(1+abs(z));
        x(i)       = z + h;
        r1         = nlp_r(x,Prob,varargin{:});
        x(i)       = x(i) + h;
        r11        = nlp_r(x,Prob,varargin{:});
        for k=ix
          d2r(i,i) = d2r(i,i)+(r(k)/h^2)*((r(k)-r1(k))+(r11(k)-r1(k)));
        end
        for j = i+1:n
           x(i)       = z + h;
           v          = x(j);
           Prob.FDVar = [i,j];
           h2         = HTol(j)*(1+abs(v));
           x(j)       = v + h2;
           r12        = nlp_r(x,Prob,varargin{:});
           x(i)       = z;
           Prob.FDVar = j;
           r2         = nlp_r(x,Prob,varargin{:});
           x(j)       = v;
           hh         = h*h2;
           % Consider only residuals for which r is nonzero
           for k=ix
             d2r(j,i) = d2r(j,i)+(r(k)/hh)*((r12(k)-r2(k))+(r(k)-r1(k)));
           end
           d2r(i,j)   = d2r(j,i);
        end
      end
   end
elseif ~isempty(Pattern) & NumDiff > 0
   % Two levels of differentiation, use r(x) values instead
   % Nonempty Prob.ConsPattern... 
   ix = find(r');
   if ~isempty(ix)
      if length(ix) == 1
         ic          = find(Pattern(ix,:));
      else
         ic          = find(any(Pattern(ix,:)));
      end
      Prob.rows      = ix;
      Prob.cols      = ic;
      for ii = 1:length(ic)
         i           = ic(ii);
         z           = x(i);
         Prob.FDVar  = i;
         h           = HTol(i)*(1+abs(z));
         x(i)        = z + h;
         r1          = nlp_r(x,Prob,varargin{:});
         x(i)        = x(i) + h;
         r11         = nlp_r(x,Prob,varargin{:});
         for k = ix
            d2r(i,i) = d2r(i,i)+(r(k)/h^2)*((r(k)-r1(k))+(r11(k)-r1(k)));
         end
         for jj = ii+1:length(ic)
            j          = ic(jj);
            x(i)       = z + h;
            v          = x(j);
            Prob.FDVar = [i,j];
            h2         = HTol(j)*(1+abs(v));
            x(j)       = v + h2;
            r12        = nlp_r(x,Prob,varargin{:});
            x(i)       = z;
            Prob.FDVar = j;
            r2         = nlp_r(x,Prob,varargin{:});
            x(j)       = v;
            hh         = h*h2;
            % Consider only residuals for which r is nonzero
            for k = ix
               d2r(j,i) = d2r(j,i)+(r(k)/hh)*((r12(k)-r2(k))+(r(k)-r1(k)));
            end
            d2r(i,j)   = d2r(j,i);
         end
      end
   end
else
   % Nonempty Prob.Jacattern... 
   ix = find(r');
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
         J_ph      = nlp_J(x,Prob,varargin{:});
         x(i)       = z;
         % Consider only residuals for which r is nonzero
         for k = ix
            d2r(ic,i) = d2r(ic,i)+ (r(k)/h)*( J_ph(k,ic)-J(k,ic) )';
         end
      end
   end
   % Make symmetric matrix
   d2r=0.5*(d2r+d2r');
end

% MODIFICATION LOG
%
% 050121 frhe File created, based on FDcHess.m
% 061117 hkh  More efficient code for safeguard against x outside simple bounds
% 090813 med  mlint check