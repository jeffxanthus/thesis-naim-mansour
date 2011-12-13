% FDHess2
%
% Numerical Hessian.
%
% The numerical approximation is computed by use of the Matlab spline
% routine (NumDiff=2) or the Spline Toolbox routines:
% csapi (NumDiff=4,SplineTol <  0), csaps (NumDiff=3,splineSmooth) or 
% spaps (NumDiff=4,SplineTol >= 0)
%
% splineTol    = Prob.optParam.splineTol
%    Should be set in the order of the noise level, default 1E-3
% splineSmooth = Prob.optParam.splineSmooth, default -1
%    0 means least squares straight line fit
%    1 means natural (variational) cubic spline interpolant
%    The transition range in [0,1] is small, suggested tries are 0.2 or 0.4
%    <0 lets the routine csaps make the choice (default)
%
% The numerical step size h is given by Prob.optParam.CentralDiff;
% 
% Sparsity pattern in Prob.HessPattern is used.
% Preprocessed efficient sparsity input in Prob.HessIx is utilized
%
% function H = FDHess2(x, Prob, gx, varargin)

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Nov 24, 1998. Last modified Aug 13, 2009.

function H = FDHess2(x, Prob, gx, varargin)

if nargin < 3
   gx = [];
end

if isempty(Prob.HessPattern)
   Pattern=[];
else
   Pattern=Prob.HessPattern;
end

x = x(:);

if isempty(gx)
   gx = nlp_g(x, Prob, varargin{:});
end

n            = Prob.N;
h            = Prob.optParam.CentralDiff;
splineSmooth = Prob.optParam.splineSmooth;
splineTol    = Prob.optParam.splineTol;
HessIx       = Prob.HessIx;
NumDiff      = abs(Prob.NumDiff);


if isempty(Pattern)
   H = zeros(n,n);
   for dim = 1:n
      z      = x(dim);
      Prob.FDVar  = dim;
      x(dim) = z + h;
      gx_ph  = nlp_g(x,Prob, varargin{:});
      x(dim) = z - h;
      gx_mh  = nlp_g(x,Prob, varargin{:});
      x(dim) = z;

      XX = [ z-h       z     z+h   ];
      for i = dim:n
         YY = [gx_mh(i) gx(i) gx_ph(i)];
      
         if NumDiff == 2       % Use spline in Matlab
            pp  = spline(XX,YY);
         elseif NumDiff == 3   % Use SPLINE Toolbox routine csaps.m
            pp  = csaps(XX,YY,splineSmooth);
         elseif splineTol >= 0 % Use SPLINE Toolbox routine spaps.m
            Bm  = spaps(XX,YY,splineTol);
            pp  = fn2fm(Bm,'pp'); % Convert B- to pp format
         else                  % Use SPLINE Toolbox routine csapi.m
            pp  = csapi(XX,YY);
         end
   
         tmp      = ppval(splineDer(pp, 1), z );
         %tmp     = fnval(fnder(pp, 1), z );
         H(i,dim) = tmp;
         H(dim,i) = tmp;
      end
   end   
   %H=0.5*(H+H');
elseif ~isempty(HessIx)
   % Assumption. Only upper rectangle is stored in Pattern
   [ix,iy,iv]=find(Pattern);
   ixS=[]; iyS=[]; ivS=[];
   mx = max(HessIx);
   for k = 1:mx
      CI         = find(HessIx==k);
      Prob.FDVar = CI;
      z          = x(CI);
      if Prob.NumDiff > 0 % Also the gradient is estimated numerically
         Prob.g_k      = gx; 
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
      x(CI)      = z + h;
      gx_ph      = nlp_g(x,Prob, varargin{:});
      x(CI)      = z - h;
      gx_mh      = nlp_g(x,Prob, varargin{:});
      x(CI )     = z;
      for j = 1:length(CI)
          iz     = find(CI(j)==iy);
          iv(iz) = (gx_ph(ix(iz))-gx(ix(iz)))/h(j);
          XX     = [ z(j)-h       z(j)     z(j)+h   ];
          for k = 1:length(iz)
              i   = ix(iz(k));
              YY  = [gx_mh(i) gx(i) gx_ph(i)];
      
              if NumDiff == 2       % Use spline in Matlab
                 pp  = spline(XX,YY);
              elseif NumDiff == 3   % Use SPLINE Toolbox routine csaps.m
                 pp  = csaps(XX,YY,splineSmooth);
              elseif splineTol >= 0 % Use SPLINE Toolbox routine spaps.m
                 Bm  = spaps(XX,YY,splineTol);
                 pp  = fn2fm(Bm,'pp'); % Convert B- to pp format
              else                  % Use SPLINE Toolbox routine csapi.m
                 pp  = csapi(XX,YY);
              end
              tmp       = ppval(splineDer(pp, 1), z );
              iv(iz(k)) = tmp;
              if i ~= CI(j)
                 ixS=[ixS;CI(j)];
                 iyS=[iyS;i];
                 ivS=[ivS;tmp];
              end
          end
      end
   end
   H=sparse([ix;ixS],[iy;iyS],[iv;ivS],n,n);
else
   % Assumption. Only upper rectangle is stored in Pattern
   [ix,iy,iv]=find(Pattern);
   ixS=[]; iyS=[]; ivS=[];
   for dim = 1:n
      % iz = find(dim==iy & dim >= ix); % Only if lower triangle is in Pattern
      iz = find(dim==iy);
      if ~isempty(iz)
         z           = x(dim);
         Prob.FDVar  = dim;
         ixz         = ix(iz);
         Prob.cols   = ixz;
         if Prob.NumDiff > 0 % Also the gradient is estimated numerically
            Prob.g_k      = gx; 
            Prob.g_k(ixz) = NaN;  % only estimate NaN elements in nlp_g
         end

         x(dim)      = z + h;
         gx_ph       = nlp_g(x,Prob, varargin{:});
         x(dim)      = z - h;
         gx_mh       = nlp_g(x,Prob, varargin{:});
         x(dim)      = z;

         XX    = [ z-h       z     z+h   ];
         for k = 1:length(iz)
            i  = ixz(k); %i=ix(iz(k));
            YY = [gx_mh(i) gx(i) gx_ph(i)];
      
            if NumDiff == 2       % Use spline in Matlab
               pp  = spline(XX,YY);
            elseif NumDiff == 3   % Use SPLINE Toolbox routine csaps.m
               pp  = csaps(XX,YY,splineSmooth);
            elseif splineTol >= 0 % Use SPLINE Toolbox routine spaps.m
               Bm  = spaps(XX,YY,splineTol);
               pp  = fn2fm(Bm,'pp'); % Convert B- to pp format
            else                  % Use SPLINE Toolbox routine csapi.m
               pp  = csapi(XX,YY);
            end
            tmp       = ppval(splineDer(pp, 1), z );
            %tmp      = fnval(fnder(pp, 1), z );
            iv(iz(k)) = tmp;
            if i ~= dim
               ixS=[ixS;dim];
               iyS=[iyS;i];
               ivS=[ivS;tmp];
            end
         end
      end
   end   
   H=sparse([ix;ixS],[iy;iyS],[iv;ivS],n,n);
end

function ppD = splineDer(pp, order)

[breaks,coeffs,L,K,D] = unmkpp(pp);
if K <= order
   ppD = mkpp([breaks(1) breaks(L+1)],zeros(D,1));
else
   Knew = K - order;
   for i=K-1:-1:Knew
       coeffs = coeffs.*repmat((i:-1:i-K+1),D*L,1);
   end
   ppD = mkpp(breaks, coeffs(:,1:Knew),D);
end

% MODIFICATION LOG
%
% 990306  hkh  Safeguard against x slightly outside bounds
% 990626  hkh  Rewritten for speed
% 000910  hkh  Written more efficiently, avoid arrays and zero multiplications
%              Use Prob.optParam.CentralDiff
% 020416  hkh  Use Prob.N for length of x, if x is longer (minimax)
% 040407  hkh  Send global index FDVar telling index of perturbed variable
% 040407  hkh  NumDiff = 2 now using standard spline, splineTol<0 => csapi
% 040407  hkh  Send Prob.FDVar telling indicies of perturbed variables
% 040407  hkh  If HessIx defined, use a more efficient method with less calls
% 040407  hkh  Set Prob.g_k with NaN, more efficient two-level differentiation
% 040407  hkh  Assume and utilize only upper triagonal of HessPattern
% 050308  hkh  spaps output must be converted to pp format using fn2fm
% 090813  med  mlint check