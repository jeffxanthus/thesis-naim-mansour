% function [g,dc] = sim_gdc(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of the gradient vector
% and the constraint gradient matrix at the same time
%
% sim_gdc calls the routine Prob.FUNCS.gdc (if nonempty)
% either as [g,dc]=feval(Prob.FUNCS.g, x) or
%           [g,dc]=feval(Prob.FUNCS.g, x, Prob,varargin)
% depending on the number of inputs
%
% The global counters n_g and n_dc are incremented

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written May 24, 2003.   Last modified Aug 14, 2006.

function [g,dc] = sim_gdc(x, Prob, varargin)

global n_g n_dc BUILDP NARG 

% Communication 
global NLP_x NLP_f NLP_pSepIndex
global NLP_c NLP_xc
global NLP_g NLP_xg NLP_dc

if ~isempty(NLP_xg) 
   if length(x)~=length(NLP_xg)
      NLP_x  = [];
      NLP_xc = [];
   elseif all(x == NLP_xg)
      g  = NLP_g;
      dc = NLP_dc;
      return
   end
end

n_g=n_g+1; n_dc=n_dc+1;

Func = Prob.FUNCS.gdc;

if isempty(Func) | (Prob.ConsDiff > 0 ) | (Prob.NumDiff > 0)
   % Here we should call a numerical difference routine
   % Check if we have computed f(x) and c(x) already (the normal case)
   fx=[];
   if ~isempty(NLP_x)
      if all(x == NLP_x) 
         fx=NLP_f;
      end
   end
   cx=[];
   if ~isempty(NLP_xc)
      if all(x == NLP_xc) 
         cx=NLP_c;
      end
   end
   BUILDP0=BUILDP;
   BUILDP=[];
   g     = [];
   dc    = [];
   [g,dc] = fdgdc(x, Prob, g, fx, dc, cx, varargin); % FD Algorithm
   NLP_x  = x;
   NLP_f  = fx;
   NLP_xc = x;
   NLP_c  = cx;
   BUILDP=BUILDP0;
else
   if isempty(NARG)
      p = xnargin(Func);
   else
      p = NARG(11);
   end
   if p > 2
      [g,dc]=feval(Func, x, Prob, varargin{:});
   elseif p > 1
      [g,dc]=feval(Func, x, Prob);
   else
      [g,dc]=feval(Func, x);
   end
   g = full(g(:));
   if Prob.CheckNaN > 0
      if any(isnan(g)) | any(isnan(dc)) 
         % There are elements set to NaN, to be estimated
         fx=[];
         if ~isempty(NLP_x)
            if all(x == NLP_x) 
               fx=NLP_f;
            end
         end
         cx=[];
         if ~isempty(NLP_xc)
            if all(x == NLP_xc) 
               cx=NLP_c;
            end
         end
         BUILDP0=BUILDP;
         BUILDP=[];
         [g,dc,fx,cx] = fdgdc(x, Prob, g, fx, dc, cx, varargin); % FD Algorithm
         NLP_x  = x;
         NLP_f  = fx;
         NLP_xc = x;
         NLP_c  = cx;
         BUILDP=BUILDP0;
      end   
   else
      g(isnan(g))   = 1E5;
      dc(isnan(dc)) = 1E5;
   end
end

NLP_xg  = x;
NLP_g   = g;
NLP_dc  = dc;

%xprint(x,'x:')
%xprint(g,'g:')
%
% Finite-Difference numerical approximation of the gradient and
% the constraint Jacobian for constraint vector c(x).
%
% function [g,dc] = fdgdc(x, Prob, g, fx, dc, cx, varargin)
% function J = FDJac(x, Prob, rxFunc, rx, varargin)
%
% fx = f(x), i.e. the current function value, if available
% cx = c(x), i.e. the current constraint vector, if available
% If g  is nonempty, estimate any elements of  g set to NaN.
% If dc is nonempty, estimate any elements of dc set to NaN.
%
% Implementation based on the algorithm FD in Practical Optimization, page 343.
%
% Sparsity pattern in Prob.ConsPattern is used.
%

function [g,dc,fx,cx] = fdgdc(x, Prob, g, fx, dc, cx, varargin)

n = Prob.N;
if isempty(g)
   ALLg = 1;
   g = zeros(n,1);
else
   ALLg = 0;
   if ~any(isnan(g)), 
      ALLg = 1;
      g = zeros(n,1);
   end
end
if isempty(Prob.ConsPattern)
   Pattern=[];
else
   Pattern=Prob.ConsPattern;
end

x = x(:);

nF = length(x);
n  = Prob.N;

x_L  = Prob.x_L(:);
if isempty(x_L)
   x_L = -inf*ones(n,1);
elseif length(x_L) < n
   x_L(length(x_L)+1:n,1) = -inf;
end
x_U  = Prob.x_U(:);
if isempty(x_U)
   x_U =  inf*ones(n,1);
elseif length(x_U) < n
   x_U(length(x_U)+1:n,1) =  inf;
end

if isempty(fx) | isempty(cx) 
   [fx,cx] = sim_fc(x, Prob, varargin{:});
end

m  = length(cx);

% gTol is a vector of relative intervals, See Pract. Opt. page 130

global gTol

if isempty(gTol)
   gTol = Prob.GradTolg;

   if isempty(gTol)
      h = Prob.optParam.DiffInt;
      if h > 0
         gTol = h*ones(n,1);
      end
   elseif gTol(1) <= 0
      gTol = [];
   else
      gTol = gTol(:);
   end
end

if isempty(gTol)
   gTol = 1E-5*ones(n,1);
end
if length(gTol) < n
   gTol = gTol(1) * ones(n,1);
end

% If close to a bound, make step in opposite direction, changing sign
ix = find(x(1:n)+gTol(1:n) > x_U | x(1:n)+gTol(1:n) < x_L);
gTol(ix)=-gTol(ix);
   
if ALLg & isempty(Pattern)
   % Compute the gradient by using the relative intervals computed for x_0.
   % All elements are estimated
   dc = zeros(m,n);
   for dim = 1:n
      Prob.FDVar    = dim;
      z             = x(dim);
      h             = gTol(dim) * (1+abs(z));
      x(dim)        = z + h;
      [fx_ph,cx_ph] = sim_fc(x, Prob, varargin{:}); % f(x+h)
      g(dim)        = (fx_ph-fx) / h;
      dc(:,dim)     = (cx_ph(:)-cx(:)) / h;
      x(dim)        = z;
   end   
elseif ALLg & ~isempty(Pattern)
   [ix,iy,iv]=find(Pattern);
   for dim = 1:n
      Prob.FDVar    = dim;
      iz = find(dim==iy);
      if ~isempty(iz)
         ixz           = ix(iz);
         Prob.rows     = ixz;
      end
      z             = x(dim);
      h             = gTol(dim) * (1+abs(z));
      x(dim)        = z + h;
      [fx_ph,cx_ph] = sim_fc(x, Prob, varargin{:}); % f(x+h)
      g(dim)        = (fx_ph-fx) / h;
      if ~isempty(iz)
         iv(iz) = (cx_ph(ixz)-cx(ixz)) / h;
      end
      x(dim)        = z;
   end   
   dc=sparse(ix,iy,iv,m,n);
elseif ~ALLg & ~isempty(Pattern)
   [ix,iy,iv]=find(Pattern);
   for dim = 1:n
      iz = find(dim==iy);
      z             = x(dim);
      h             = gTol(dim) * (1+abs(z));
      x(dim)        = z + h;
      if isnan(g(dim)) | ~isempty(iz)
         Prob.FDVar    = dim;
         if ~isempty(iz)
            ixz           = ix(iz);
            Prob.rows     = ixz;
         end
         [fx_ph,cx_ph] = sim_fc(x, Prob, varargin{:}); % f(x+h)
      end
      if isnan(g(dim))
         g(dim)        = (fx_ph-fx) / h;
      end
      if ~isempty(iz)
         iv(iz) = (cx_ph(ixz)-cx(ixz)) / h;
      end
      x(dim)        = z;
   end   
   dc=sparse(ix,iy,iv,m,n);
elseif ~ALLg & isempty(Pattern)
   dc = zeros(m,n);
   for dim = 1:n
      Prob.FDVar    = dim;
      z             = x(dim);
      h             = gTol(dim) * (1+abs(z));
      x(dim)        = z + h;
      [fx_ph,cx_ph] = sim_fc(x, Prob, varargin{:}); % f(x+h)
      if isnan(g(dim))
         g(dim)        = (fx_ph-fx) / h;
      end
      dc(:,dim)     = (cx_ph-cx) / h;
      x(dim)        = z;
   end   
end

% MODIFICATION LOG
%
% 030524 hkh  Written
% 040304 hkh  Safe guard for NaN
% 040318 hkh  Return fx,cx from fdgdc, otherwise fx,cx might be set []
% 040411 hkh  Send Prob.FDVar telling indicies of perturbed variables
% 040411 hkh  Prob.rows-constraints needed to be computed, if Pattern not []
% 040413 hkh  Add test to avoid going outside simple bounds
% 060814 med  FUNCS used for callbacks instead