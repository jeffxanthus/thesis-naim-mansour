% FDJac2
%
% Compute an approximation of the Jacobian of the residuals in nonlinear
% least square problem or the constraint Jacobian for constraint vector c(x).
% The numerical approximation is computed by use of the Matlab spline
% routine (NumDiff=2) or the Spline Toolbox routines
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
% function [J,rx] = FDJac2(x, Prob, rxFunc, rx, varargin)
%
% rx = r(x) or c(x), i.e. the current residual or constraint vector,
% if available, otherwise left empty.
%
% Sparsity pattern in Prob.JacPattern or Prob.ConsPattern is used.
% Preprocessed efficient sparsity input in Prob.ConIx is utilized

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Oct 27, 1998. Last modified Aug 13, 2009.

function [J,rx] = FDJac2(x, Prob, rxFunc, rx, varargin)

if nargin < 4
   rx = [];
   if nargin < 3
      rxFunc = 'nlp_r';
   end
end

if strcmp(rxFunc,'nlp_r')
   NLP = 1;
elseif strcmp(rxFunc,'nlp_c')
   NLP = 2;
else
   NLP = 0;
end
if NLP==1 | NLP==0
   NumDiff = abs(Prob.NumDiff);
   if isempty(Prob.JacPattern)
      Pattern = [];
      ConIx = [];
   else
      Pattern = Prob.JacPattern;
      ConIx = Prob.JacIx;
   end
elseif NLP==2
   NumDiff = abs(Prob.ConsDiff);
   if isempty(Prob.ConsPattern)
      Pattern = [];
      ConIx = [];
   else
      Pattern = Prob.ConsPattern;
      ConIx = Prob.ConIx;
   end
end

x = x(:);
splineSmooth = Prob.optParam.splineSmooth;
splineTol    = Prob.optParam.splineTol;

if isempty(rx)
   if NLP==1
      rx = nlp_r(x, Prob, varargin{:});
   elseif NLP==2
      rx = nlp_c(x, Prob, varargin{:});
   else
      rx = feval(rxFunc, x, Prob, varargin{:});
   end
end

n = Prob.N;
m = length(rx);
h = Prob.optParam.CentralDiff;

if isempty(Pattern)
   if Prob.LargeScale
      J = sparse(m,n);
   else
      J = zeros(m,n);
   end

   for dim = 1:n
      z = x(dim);
   
      x(dim) = z + h;
      Prob.FDVar = dim;
      Prob.cols  = dim;
      if strcmp(rxFunc,'nlp_r')
         rx_ph   = nlp_r(x,Prob, varargin{:});
         x(dim)  = z - h;
         rx_mh   = nlp_r(x,Prob, varargin{:});
      elseif strcmp(rxFunc,'nlp_c')
         rx_ph   = nlp_c(x,Prob, varargin{:});
         x(dim)  = z - h;
         rx_mh   = nlp_c(x,Prob, varargin{:});
      else
         rx_ph   = feval(rxFunc,x,Prob, varargin{:});
         x(dim)  = z - h;
         rx_mh   = feval(rxFunc,x,Prob, varargin{:});
      end
      x(dim) = z;
      
      XX = [ z-h       z     z+h   ];
      for i = 1:m
         YY = [rx_mh(i) rx(i) rx_ph(i)];
         
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
          
         J(i,dim) = ppval(splineDer(pp, 1), z );
         %J(i,dim) = fnval(fnder(pp, 1), z );
      end
   end   
elseif ~isempty(ConIx)
   [ix,iy,iv]=find(Pattern);
   iv = double(iv);
   mx = max(ConIx);
   for k = 1:mx
      CI         = find(ConIx==k);
      nCI        = length(CI);
      Prob.FDVar = CI;
      Prob.cols  = CI;
      if length(CI) == 1
         Prob.rows  = ix(find(CI==iy));
      else
         Prob.rows  = ix(find(any([ones(length(iy),1)*CI==iy*ones(1,nCI)]')));
      end
      z          = x(CI);
      x(CI)      = z + h;
      if NLP==1
         rx_ph   = nlp_r(x,Prob, varargin{:});
         x(CI)   = z - h;
         rx_mh   = nlp_r(x,Prob, varargin{:});
      elseif NLP==2
         rx_ph   = nlp_c(x,Prob, varargin{:});
         x(CI)   = z - h;
         rx_mh   = nlp_c(x,Prob, varargin{:});
      else
         rx_ph   = feval(rxFunc,x,Prob, varargin{:});
         x(CI)   = z - h;
         rx_mh   = feval(rxFunc,x,Prob, varargin{:});
      end
      x(CI)      = z;
      for j = 1:length(CI)
          iz = find(CI(j)==iy);
          XX = [ z(j)-h       z(j)     z(j)+h   ];
          for kk = 1:length(iz)
             i=ix(iz(kk));
             YY = [rx_mh(i) rx(i) rx_ph(i)];
         
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
             iv(iz(kk)) = ppval(splineDer(pp, 1), z(j) );
             %iv(iz(kk)) = fnval(fnder(pp, 1), z(j) );
          end
      end
   end
   J=sparse(ix,iy,iv,m,n);
else
   [ix,iy,iv]=find(Pattern);
   iv = double(iv);

   for dim = 1:n

       iz         = find(dim==iy);
       ixz        = ix(iz);
       Prob.rows  = ixz;
       if ~isempty(iz)
          z          = x(dim);
          Prob.FDVar = dim;
          Prob.cols  = dim;
          x(dim)     = z + h;
          if NLP==1
             rx_ph   = nlp_r(x,Prob, varargin{:});
             x(dim)  = z - h;
             rx_mh   = nlp_r(x,Prob, varargin{:});
          elseif NLP==2
             rx_ph   = nlp_c(x,Prob, varargin{:});
             x(dim)  = z - h;
             rx_mh   = nlp_c(x,Prob, varargin{:});
          else
             rx_ph   = feval(rxFunc,x,Prob, varargin{:});
             x(dim)  = z - h;
             rx_mh   = feval(rxFunc,x,Prob, varargin{:});
          end
          x(dim)     = z;

          XX         = [ z-h       z     z+h   ];
          for k = 1:length(iz)
              i      = ixz(k);
              YY     = [rx_mh(i) rx(i) rx_ph(i)];
         
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
          
              iv(iz(k)) = ppval(splineDer(pp, 1), z );
              %iv(iz(k)) = fnval(fnder(pp, 1), z );
          end
       end
   end   
   J=sparse(ix,iy,iv,m,n);
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
% 981027  hkh  Return empty J in case of bound error on x. Spelling errors
% 981028  hkh  Rewriting fdng2 as general routine for both c and r.
%              Add rxFunc as input. 
% 981110  mbk  Added call to SPLINE Toolbox routines csaps.m and spaps.m
%              if Prob.NumDiff equals 3 or 4.
% 981124  mbk  Use of feval when calling SPLINE Toolbox routines.
% 990306  hkh  Safeguard against x slightly outside bounds
% 000910  hkh  Written more efficiently, avoid arrays and zero multiplications
%              Use Prob.optParam.CentralDiff
% 020416  hkh  Use Prob.N for length of x, if x is longer (minimax)
% 031214  hkh  Use J = sparse(m.n) if Prob.LargeScale = 1
% 040125  hkh  Set double(iv) to avoid treatment at logical array in Matlab 6.5
% 040331  hkh  Return rx as second parameter
% 040406  hkh  If ConIx defined, use a more efficient method with less calls
% 040406  hkh  Send Prob.FDVar telling indicies of perturbed variables
% 040406  hkh  NumDiff = 2 now using standard spline, splineTol<0 => csapi
% 040410  hkh  Send Prob.cols and Prob.rows
% 040505  hkh  Incorrect use of mx as length of CI
% 050124  hkh  ConIx=Prob.JacIx must be set for NLLS problems
% 050308  hkh  spaps output must be converted to pp format using fn2fm
% 090813  med  mlint check