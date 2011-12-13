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
      parfor i = 1:n
        xP           = x;
        z            = xP(i);
        %Prob=       = Prob;
        %ProbP.FDVar = i;
        h            = HTol(i)*(1+abs(z));
        xP(i)        = z + h;
        dc_ph        = nlp_dc(xP,Prob,varargin{:});
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
      %Prob.rows    = ix;
      c            = nlp_c(x,Prob,varargin{:});
      parfor i = 1:n
         xP            = x;
         z             = xP(i);
         % Change Prob to ProbP in 4 calls to nlp_c if using the next 2 lines
         %ProbP        = Prob;
         %Prob.FDVar   = i;
         h             = HTol(i)*(1+abs(z));
         xP(i)         = z + h;
         c1            = nlp_c(xP,Prob,varargin{:});
         xP(i)         = xP(i) + h;
         c11           = nlp_c(xP,Prob,varargin{:});
         d2cP(i).di    = sum((lam(ix)/h^2).*((c(ix)-c1(ix))+(c11(ix)-c1(ix))));
	 d             = zeros(1,n-i)
         for j = i+1:n
            xP(i)        = z + h;
            v            = xP(j);
            %Prob.FDVar  = [i,j];
            h2           = HTol(j)*(1+abs(v));
            xP(j)        = v + h2;
            c12          = nlp_c(xP,Prob,varargin{:});
            xP(i)        = z;
            %Prob.FDVar  = j;
            c2           = nlp_c(xP,Prob,varargin{:});
            xP(j)        = v;
            hh           = h*h2;
            % Consider only constraints for which lam is nonzero
            d(j-i)       = sum((lam(ix)/hh).*((c12(ix)-c2(ix))+(c(ix)-c1(ix))));
          end
          d2cP(i).d         = d;
      end
      for i = 1:n
          d2c(i,i)      = d2cP(i).di;
          d             = d2cP(i).d;
          d2c(i,i+1:n)  = d;
          d2c(i+1:n,i)  = d';
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
      HessPatt       = sparse(n,n);
      nic            = length(ic);
      % Prob.rows      = ix;
      c              = nlp_c(x,Prob,varargin{:});
      % Prob.cols      = ic;
      parfor ii = 1:nic
         xP           = x;
         i            = ic(ii);
         z            = xP(i);
         % Change Prob to ProbP in 4 calls to nlp_c if using the ProbP 3 lines
         %ProbP.FDVar = i;
         h            = HTol(i)*(1+abs(z));
         xP(i)        = z + h;
         c1           = nlp_c(xP,Prob,varargin{:});
         xP(i)        = xP(i) + h;
         c11          = nlp_c(xP,Prob,varargin{:});
         d2cP(ii).di  = sum((lam(ix)/h^2).*((c(ix)-c1(ix))+(c11(ix)-c1(ix))));
	 d            = zeros(1,nic-ii);
         for jj = ii+1:nic
            j            = ic(jj);
            xP(i)        = z + h;
            v            = xP(j);
            %ProbP.FDVar = [i,j];
            h2           = HTol(j)*(1+abs(v));
            xP(j)        = v + h2;
            c12          = nlp_c(xP,Prob,varargin{:});
            xP(i)        = z;
            %ProbP.FDVar = j;
            c2           = nlp_c(xP,Prob,varargin{:});
            xP(j)        = v;
            hh           = h*h2;
            % Consider only constraints for which lam is nonzero
            d(jj-ii)     = sum((lam(ix)/hh).*((c12(ix)-c2(ix))+(c(ix)-c1(ix))));
            %for k=ix
            %    d(jj-ii) = d(jj-ii)+(lam(k)/hh)*((c12(k)-c2(k))+(c(k)-c1(k)));
            %end
         end
         d2cP(ii).d      = d;
      end
      for ii = 1:nic
          i             = ic(ii);
          d2c(i,i)      = d2cP(ii).di;
          d             = d2cP(ii).d;
	  idx           = ic(ii+1:nic);
          d2c(i,idx)    = d;
          d2c(idx,i)    = d';
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
      nic           = length(ic);
      d2cP           = sparse(nic,nic); % Compressed matrix
      %Prob.rows     = ix;
      %Prob.cols     = ic;
      parfor ii = 1:nic
          i            = ic(ii);
          xP           = x;
          z            = xP(i);
          % Change Prob to ProbP in call to nlp_dc if using ProbP 2 lines, Prob 2 lines
          %ProbP       = Prob;
          %ProbP.FDVar = i;
          z            = xP(i);
          h            = HTol(i)*(1+abs(z));
          xP(i)        = z + h;
          dc_ph        = nlp_dc(xP,Prob,varargin{:});
          % Consider only constraints for which lam is nonzero
          for k = ix
             d2cP(:,ii) = d2cP(:,ii)+ (lam(k)/h)*( dc_ph(k,ic)-dc(k,ic) )';
          end
      end
   end
   % Make symmetric matrix
   d2cP          = 0.5*(d2cP+d2cP');
   % Expand to full symmetric matrix
   for ii = 1:nic
       i         = ic(ii);
       d2c(ic,i) = d2cP(:,ii);
   end
end

% MODIFICATION LOG
%
% 110727 hkh  Written from FDcHess
% 110727 hkh  Use 5*HTol for 2 levels of differentiation
