% FDHess_par
%
% Numerical approximation of the Hessian matrix
% This version is using parfor and the Parallel Computing Toolbox
%
% Implementation based on the algorithm FD in Practical Optimization, page 343.
%
% function H = FDHess_par(x, Prob, gx, varargin)
%
% Sparsity pattern in Prob.HessPattern is used.
% Preprocessed efficient sparsity input in Prob.HessIx is utilized

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Nov 24, 1998.  Last modified July 27, 2011.

function H = FDHess_par(x, Prob, gx, varargin)

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
   % Case 1: NumDiff < 0, i.e. one level of differentiation, no pattern defined 
   parfor dim = 1:n
      xP           = x;
      % Change Prob to ProbP in nlp_g call if using the 2 next lines
      %ProbP       = Prob;
      %ProbP.FDVar = dim;
      z            = xP(dim);
      h            = HTol(dim)*(1+abs(z));
      xP(dim)      = z + h;
      gx_ph        = nlp_g(xP,Prob, varargin{:});
      H(:,dim)     = (gx_ph-gx)/h;
   end   
   H=0.5*(H+H'); % Make perfectly symmetric
elseif isempty(Pattern)
   % Case 2: NumDiff > 0, i.e. two levels of differentiation, 
   % use f(x) values instead. No pattern defined
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
   parfor i = 1:n
      w            = zeros(1,n);
      xP           = x;
      z            = xP(i);
      % Change Prob to ProbP in 4 calls to nlp_f if using the 2 next lines and 2 other ProbP
      %ProbP       = Prob;
      %ProbP.FDVar = i;
      h            = HTol(i)*(1+abs(z));
      xP(i)        = z + h;
      f1           = nlp_f(xP,Prob,varargin{:});
      xP(i)        = xP(i) + h;
      f11          = nlp_f(xP,Prob,varargin{:});
      w(i)         = ((f-f1)+(f11-f1))/h^2;
      for j = i+1:n
         xP(i)        = z + h;
         v            = xP(j);
         %ProbP.FDVar = [i,j];
         h2           = HTol(j)*(1+abs(v));
         xP(j)        = v + h2;
         f12          = nlp_f(xP,Prob,varargin{:});
         xP(i)        = z;
         %ProbP.FDVar = j;
         f2           = nlp_f(xP,Prob,varargin{:});
         xP(j)        = v;
         w(j)         = ((f12-f2)+(f-f1))/(h*h2);
      end
      HP(i).w = w;
   end
   for i = 1:n
       z        = HP(i).w(i:n);
       z
       H(i,i:n) = z;
       H(i:n,i) = z';
   end
elseif ~isempty(HessIx)
   % Case 3: HessIx defined (i.e. also Pattern)
   % Two levels of differentiation, use f(x) values instead
   % Assumption. Only upper rectangle is stored in Pattern
   ixS=[]; iyS=[]; ivS=[];
   mx = max(HessIx);
   parfor k = 1:mx
      xP            = x;
      CI            = find(HessIx==k);
      % Change Prob to ProbP in call to nlp_g if using the 2 next lines and 4 other ProbP
      ProbP         = Prob;
      ProbP.FDVar   = CI;
      z             = xP(CI);
      h             = HTol(CI).*(1+abs(z));
      if NumDiff > 0 % Also the gradient is estimated numerically
         ProbP.g_k  = gx; 
         for j = 1:length(CI)
             % only estimate NaN elements in nlp_g
             ProbP.g_k(ix(find(CI(j)==iy))) = NaN;  
         end
         ProbP.cols = find(isnan(ProbP.g_k));
      else
         if length(CI) == 1
            ProbP.cols = ix(find(CI==iy));
         else
            ProbP.cols = ix(find(any([ones(length(iy),1)*CI==iy*ones(1,mx)]')));
         end
      end
      xP(CI)        = z + h;
      gx_ph         = nlp_g(xP,ProbP, varargin{:});
      HP(k).g       = gx_ph;
   end
   for k = 1:mx
      CI            = find(HessIx==k);
      z             = x(CI);
      h             = HTol(CI).*(1+abs(z));
      gx_ph         = HP(k).g;
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
   % Case 4: NumDiff > 0, nonempty Pattern
   % Two levels of differentiation, use f(x) values instead
   % Assumption. Only upper rectangle is stored in Pattern
   [ix,iy,iv]=find(Pattern);
   % Increase tolerance by factor 10
   HTol            = 10*HTol;
   f               = nlp_f(x,Prob,varargin{:});
   iD              = unique(iy);
   niD             = length(iD);
   parfor ii = 1:niD
      xP           = x;
      % Change Prob to ProbP in 4 calls to nlp_f if using the 2 next lines and 2 other ProbP
      %ProbP       = Prob;
      %ProbP.FDVar = i;
      i            = iD(ii);
      iz           = find(i==iy);     % Which elements have 2nd index i
      ixz          = ix(iz);          % Which row elements that have 2nd index i
      m            = length(ixz)-1;
      z            = xP(i);
      h            = HTol(i)*(1+abs(z));
      xP(i)        = z + h;
      f1           = nlp_f(xP,Prob,varargin{:});
      xP(i)        = xP(i) + h;
      f11          = nlp_f(xP,Prob,varargin{:});
      y            = zeros(m+1,1);
      y(end)       = ((f-f1)+(f11-f1))/h^2;
      for jj = 1:m
         j            = ixz(jj);
         xP(i)        = z + h;
         v            = xP(j);
         %ProbP.FDVar = [i,j];
         h2           = HTol(j)*(1+abs(v));
         xP(j)        = v + h2;
         f12          = nlp_f(xP,Prob,varargin{:});
         xP(i)        = z;
         %ProbP.FDVar = j;
         f2           = nlp_f(xP,Prob,varargin{:});
         xP(j)        = v;
         y(jj)        = ((f12-f2)+(f-f1))/(h*h2);
      end
      xP(i)           = z;
      HP(ii).iz       = iz;
      HP(ii).y        = y;
      HP(ii).ixz      = ixz;
   end
   % The following method is faster than the one in comments
   ixS                = zeros(niD,1); 
   iyS                = zeros(niD,1); 
   ivS                = zeros(niD,1);
   nS                 = 0;
   for ii = 1:niD
       i              = iD(ii);
       iz             = HP(ii).iz;
       y              = HP(ii).y;
       ixz            = HP(ii).ixz;
       iv(iz)         = y;
       m              = length(ixz)-1;
       ixS(nS+1:nS+m) = i;
       iyS(nS+1:nS+m) = ixz(1:m);
       ivS(nS+1:nS+m) = y(1:m);
       nS             = nS + m;
   end
   H=sparse([ix;ixS(1:nS)],[iy;iyS(1:nS)],[iv;ivS(1:nS)],n,n);
   %H                        = sparse(n,n);
   %for ii = 1:niD
   %    i                    = iD(ii);
   %    ixz                  = HP(ii).ixz;
   %    y                    = HP(ii).y;
   %    H(ixz,i)             = y;
   %    if length(y) > 1
   %       H(i,ixz(1:end-1)) = y(1:end-1)';
   %    end
   %end
else
   % Case 5: NumDiff < 0, only 1 level of differentiation, nonempty Pattern
   % Assumption. Only upper rectangle is stored in Pattern
   [ix,iy,iv]        = find(Pattern);
   iD                = unique(iy);
   niD               = length(iD);
   H                 = sparse(n,n);
   HP                = sparse(n,niD);
   parfor i = 1:niD
      xP             = x;       
      dim            = iD(i);
      % Change Prob to ProbP in 1 call to nlp_g if using the next 2 lines
      %ProbP        = Prob;
      %ProbP.FDVar    = dim;
      z              = xP(dim);
      h              = HTol(dim)*(1+abs(z));
      xP(dim)        = z + h;
      gx_ph          = nlp_g(xP,Prob, varargin{:});
      HP(:,i)        = sparse((gx_ph-gx)/h);
   end   
   % The following is slower than the next
   %for i = 1:niD
   %    dim             = iD(i);
   %    H(dim:n,dim)    = HP(dim:n,i);
   %    H(dim,dim+1:n)  = H(dim+1:n,dim)';
   %    %H(dim,dim:n)    = H(dim:n,dim)';
   %end
   for i = 1:niD
       dim           = iD(i);
       z             = HP(dim:n,i);
       H(dim:n,dim)  = z;
       H(dim,dim:n)  = z';
   end
end

% MODIFICATION LOG
%
% 110725  hkh  Made parfor version out of FDHess
% 110725  hkh  Save time avoiding sending FDVar info to nlp_f/nlp_g
% 110727  hkh  Use 10*HTol is 2 levels of differentiation

