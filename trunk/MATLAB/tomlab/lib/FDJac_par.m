% FDJac_par
%
% Numerical approximation of the Jacobian of the residuals in nonlinear
% least squares problem or the constraint Jacobian for constraint vector c(x).
% This version is using parfor and the Parallel Computing Toolbox
%
% function [J,rx] = FDJac_par(x, Prob, rxFunc, rx, varargin)
%
% rx = r(x) or c(x), i.e. the current residual or constraint vector,
% if available
%
% Implementation based on the algorithm FD in Practical Optimization, page 343.
%
% Sparsity pattern in Prob.JacPattern or Prob.ConsPattern is used.
% Preprocessed efficient sparsity input in Prob.ConIx is utilized

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Nov 23, 1998.  Last modified July 25, 2011.

function [J,rx] = FDJac_par(x, Prob, rxFunc, rx, varargin)

if nargin < 4
   rx = [];
   if nargin < 3
      rxFunc = 'nlp_r';
   end
end

if strcmp(rxFunc,'nlp_r')
   NLP=1;
elseif strcmp(rxFunc,'nlp_c')
   NLP=2;
else
   NLP=0;
end
if NLP==1 | NLP==0
   if isempty(Prob.JacPattern)
      Pattern = [];
      ConIx   = [];
   else
      Pattern = Prob.JacPattern;
      ConIx   = Prob.JacIx;
   end
elseif NLP==2
   if isempty(Prob.ConsPattern)
      Pattern = [];
      ConIx   = [];
   else
      Pattern = Prob.ConsPattern;
      ConIx   = Prob.ConIx;
   end
end

x = x(:);
n  = Prob.N;

if isempty(rx)
   rx = rEval(rxFunc, x, Prob, NLP, varargin{:});
end
if isempty(rx)
   % If there are no nonlinear constraints then dc should be [].
   J = []; return;
end

fx = 0.5*(rx'*rx);
m  = length(rx);

% JTol is vector of relative intervals, See Pract. Opt. page 130

global JTol

if isempty(JTol)
   JTol = Prob.GradTolJ;

   if isempty(JTol)

      h = Prob.optParam.DiffInt;
      if h > 0
         JTol = h*ones(n,1);
      end
   elseif JTol(1) <= 0
      JTol = [];
   else
      JTol=JTol(:);
   end
end

if isempty(JTol)

   % First time this routine is called in the current optimization.
   % Run the FD algorithm to determine the relative intervals.
   %
   % This procedure is inefficient for sparse J, and does not take into
   % account the pattern in JacPattern

   % rho = 10;
   eps_R = eps;
   eps_A = eps_R*(1 + abs(fx));
   eps_phi = 1E-14; % zero toleranse for phi_F, phi_B and PHI.
   M = 1E10; % Large number.
   K = 6;
   % eta = 1, omega = 1

   PriLev = Prob.PriLevOpt;

   JTol = zeros(n,1);
   if Prob.LargeScale
      J = sparse(m,n);
   else
      J = zeros(m,n);
   end

   for dim = 1:n
      z          = x(dim);
      Prob.FDVar = dim;
      % Flags:
      FD3 = 1; FD4 = 1; FD5 = 1; FD6 = 1;

      % FD1: Initialization
      % h_line  = 2*(eta + abs(z))*sqrt(eps_A/(omega+abs(fx)));
      h_line  = 2*(1 + abs(z))*sqrt(eps_A/(1+abs(fx)));
      h(1)    = 10*h_line; % h(i+1) corresponds to h(i) in Pract. Opt.
      k       = 0;

      x(dim)  = z + h(1);
      rx_ph   = rEval(rxFunc,x,Prob,NLP,varargin{:}); % r(x+h)
      fx_ph   = 0.5*(rx_ph'*rx_ph);
      x(dim)  = z - h(1);
      rx_mh   = rEval(rxFunc,x,Prob,NLP,varargin{:}); % f(x-h)
      x(dim)  = z;
      fx_mh   = 0.5*(rx_mh'*rx_mh);
      phi_F   = (fx_ph-fx) / h(1);
      phi_B   = (fx-fx_mh) / h(1);
      phi_C   = (fx_ph-fx_mh) / (2*h(1));
      PHI     = (fx_ph - 2*fx + fx_mh) / ( h(1)^2 );
      if abs(phi_F) <= eps_phi
         C_phi_F = M;
      else
         C_phi_F = 2*eps_A / ( h(1) * abs(phi_F) );
      end
      if abs(phi_B) <= eps_phi
         C_phi_B = M;
      else
         C_phi_B = 2*eps_A / ( h(1) * abs(phi_B) );
      end
      if abs(PHI) <= eps_phi
         C_PHI = M;
      else
         C_PHI   = 4*eps_A / ( h(1)^2 * abs(PHI) );
      end
      h_s = -1;


      % FD2: Decide whether to accept the initial interval
      if max(C_phi_F,C_phi_B) <= 0.1
         h_s = h(1);
      end
      if ( C_PHI >= 0.001 ) & ( C_PHI <= 0.1 )
         h_phi = h(1);
         FD3 = 0; FD4 = 0;  % GOTO STEP FD5
      elseif C_PHI < 0.001
         FD3 = 0;           % GOTO STEP FD4
      end

      % FD3: Increase h
      while FD3 & (k < K)
         k = k + 1;
         h(k+1) = 10*h(k);
         % Compute the associated finite-difference estimates and their
         % relative errors:
         x(dim)  = z + h(k+1);
         rx_ph   = rEval(rxFunc,x,Prob,NLP,varargin{:});
         fx_ph   = 0.5*(rx_ph'*rx_ph);
         x(dim)  = z - h(k+1);
         rx_mh   = rEval(rxFunc,x,Prob,NLP,varargin{:});
         x(dim)  = z;
         fx_mh   = 0.5*(rx_mh'*rx_mh);
         phi_F   = (fx_ph-fx) / h(k+1);
         phi_B   = (fx-fx_mh) / h(k+1);
         PHI     = (fx_ph - 2*fx + fx_mh) /  h(k+1)^2 ;
         if abs(phi_F) <= eps_phi
            C_phi_F = M;
         else
            C_phi_F = 2*eps_A / ( h(k+1) * abs(phi_F) );
         end
         if abs(phi_B) <= eps_phi
            C_phi_B = M;
         else
            C_phi_B = 2*eps_A / ( h(k+1) * abs(phi_B) );
         end
         if abs(PHI) <= eps_phi
            C_PHI = M;
         else
            C_PHI   = 4*eps_A / ( h(k+1)^2 * abs(PHI) );
         end
         if (h_s < 0) & (max(C_phi_F,C_phi_B) <= 0.1)
            h_s = h(k+1);
         end
         if C_PHI <= 0.1
            h_phi = h(k+1);
            FD3 = 0; FD4 = 0; % GOTO STEP FD5
         end
      end
      if k == K
         FD4 = 0; FD5 = 0; % GOTO STEP FD6
      end

      while FD4 & (k < K)
         % FD4: Decrease h
         k = k + 1;
         h(k+1) = h(k)/10;
         % Compute the associated finite-difference estimates and their
         % relative errors:
         x(dim)  = z + h(k+1);
         rx_ph   = rEval(rxFunc,x,Prob,NLP,varargin{:});
         fx_ph   = 0.5*(rx_ph'*rx_ph);
         x(dim)  = z - h(k+1);
         rx_mh   = rEval(rxFunc,x,Prob,NLP,varargin{:});
         x(dim)  = z;
         fx_mh   = 0.5*(rx_mh'*rx_mh);
         phi_F   = (fx_ph-fx) / h(k+1);
         phi_B   = (fx-fx_mh) / h(k+1);
         PHI     = (fx_ph - 2*fx + fx_mh) / ( h(k+1)^2);
         if abs(phi_F) <= eps_phi
            C_phi_F = M;
         else
            C_phi_F = 2*eps_A / ( h(k+1) * abs(phi_F) );
         end
         if abs(phi_B) <= eps_phi
            C_phi_B = M;
         else
            C_phi_B = 2*eps_A / ( h(k+1) * abs(phi_B) );
         end
         if abs(PHI) <= eps_phi
            C_PHI = M;
         else
            C_PHI   = 4*eps_A / ( h(k+1)^2 * abs(PHI) );
         end
         if C_PHI > 0.1
            h_phi = h(k);
            FD4 = 0; % GOTO STEP FD5
         else
            if max(C_phi_F,C_phi_B) < 0.1
               h_s = h(k+1);
            end
            if ( C_PHI >= 0.001 ) & ( C_PHI <= 0.1 )
               h_phi = h(k+1);
               FD4 = 0; % GOTO STEP FD5
            end
         end
      end
      if k == K
         FD5 = 0; % GOTO STEP FD6
      end
      % HKH safe guard
      if PHI == 0, PHI = 1E-7; end

      if FD5
         % FD5: Compute the estimate of the optimal interval
         h_F = 2 * sqrt( eps_A / abs(PHI) );
         x(dim)  = z + h_F;
         rx_ph   = rEval(rxFunc,x,Prob,NLP,varargin{:});
         x(dim)  = z;
         rx_phh  = rx_ph;
         fx_ph   = 0.5*(rx_ph'*rx_ph);
         phi     = (fx_ph-fx) / h_F;

         % Compute the estimated error bound
         E_F     = h_F*abs(PHI)/2 + 2*eps_A/h_F;

         % Compute the difference between phi and phi_C(h_phi)
         x(dim)  = z + h_phi;
         rx_ph   = rEval(rxFunc,x,Prob,NLP,varargin{:}); % f(x+h)
         fx_ph   = 0.5*(rx_ph'*rx_ph);
         x(dim)  = z - h_phi;
         rx_mh   = rEval(rxFunc,x,Prob,NLP,varargin{:}); % f(x-h)
         x(dim)  = z;
         fx_mh   = 0.5*(rx_mh'*rx_mh);
         phi_C   = (fx_ph-fx_mh) / (2*h_phi);
         E_line  = abs(phi - phi_C);

         if max(E_F,E_line) <= 0.5*abs(phi)
            FD6 = 0; % the algorithm terminates successfully
         else
            FD6 = 0;
            if PriLev > 0
               disp('Termination with an error condition in FD5 in FDJac.m')
            end
         end
      end

      if FD6
         % FD6: Check unsatisfactory cases
         if PriLev > 2
            disp(' FD6: Check cases when no satisfactory interval was found');
         end
         if h_s < 0
            if PriLev > 2
               disp(' f appear to be nearly constant')
            end
            h_F = h_line;
            rx_phh = rx; % This will give J(:,dim)=0;
            %  phi = 0;
            %  PHI = 0;
            %  E_F = 0;
         elseif ( C_PHI > 0.1 ) & ( h_s > 0 )
            if PriLev > 2
               disp(' f appears to be odd or nearly linear')
            end
            h_F    = h_s;
            x(dim) = z + h_F;
            rx_phh = rEval(rxFunc,x,Prob,NLP,varargin{:});
            x(dim) = z;
            %  fx_ph   = 0.5*rx_ph'*rx_ph;
            %  phi   = (fx_ph-fx) / h_F;
            %  PHI   = 0;
            %  E_F   = 2*eps_A/h_F;
         else
            if PriLev > 2
               disp(' d2f/dx2 appears to be increasing rapidly as h decreases')
            end
            h_F    = h(K+1);
            x(dim) = z + h_F;
            rx_phh = rEval(rxFunc,x,Prob,NLP,varargin{:});
            x(dim) = z;
            %  fx_ph   = 0.5*rx_ph'*rx_ph;
            %  phi   = (fx_ph-fx) / h_F;
            %  rx_mh = rEval(rxFunc,x-h_F*E,Prob,NLP,varargin{:}); % f(x-h)
            %  fx_mh   = 0.5*rx_mh'*rx_mh;
            %  PHI   = (fx_ph - 2*fx + fx_mh) / ( h_F^2 );
            %  E_F   = h_F*abs(PHI)/2 + 2*eps_A/h_F;
         end
      end
      J(:,dim) = (rx_phh-rx)./h_F;
      JTol(dim) = h_F/(1+abs(z));
   end % for
else

   % Do not run the FD algorithm, compute the gradient by using the
   % relative intervals computed for x_0.
   if length(JTol) < n
      JTol(1:n) = JTol(1);
      JTol      = JTol(:);
   end

   % Safeguard against x outside simple bounds
   ixx = zeros(n,1);
   ixx(1:n)=ixx(1:n)+(x(1:n)+JTol(1:n).*(1+abs(x(1:n))) < Prob.x_L(1:n));
   ixx(1:n)=ixx(1:n)+(x(1:n)+JTol(1:n).*(1+abs(x(1:n))) > Prob.x_U(1:n));
   % If violating bound, make step in opposite direction, changing sign
   ixx = ixx(ixx>0);
   JTol(ixx)=-JTol(ixx);

   if ~isempty(ConIx)
      % Case 1. ConIX defined prob pattern, ConsDiff > 10
      [ix,iy,iv]    = find(Pattern);
      iv            = double(iv);
      mx            = max(ConIx);
      for k = 1:mx
         xP             = x;
         CI             = find(ConIx==k);
         nCI            = length(CI);
         % Change Prob to ProbP in rxFunc call if using the 3 next lines
         %ProbP         = Prob;
         %ProbP.FDVar   = CI;
         %ProbP.cols    = CI;
         %if nCI == 1
         %   ProbP.rows = ix(find(CI==iy));
         %else
         %   %Prob.rows = ix(find(any([ones(length(iy),1)*CI==iy*ones(1,mx)]')));
         %   ProbP.rows = ix(find(any([ones(length(iy),1)*CI==iy*ones(1,nCI)]')));
         %end
         z              = xP(CI);
         h              = JTol(CI).*(1+abs(z));
         xP(CI)         = z + h;
         rx_ph          = rEval(rxFunc,xP,Prob,NLP,varargin{:});
	 HP(k).rx       = rx_ph; 
      end
      for k = 1:mx
         x(CI )         = z;
         CI             = find(ConIx==k);
         nCI            = length(CI);
         z              = xP(CI);
         h              = JTol(CI).*(1+abs(z));
	 rx_ph          = HP(k).rx;
         for j = 1:nCI
             iz     = find(CI(j)==iy);
             iv(iz) = (rx_ph(ix(iz))-rx(ix(iz)))/h(j);
         end
      end
      J = sparse(ix,iy,iv,m,n);
   elseif isempty(Pattern)
      % Case 2. No pattern defined
      if Prob.LargeScale
         J = sparse(m,n);
      else
         J = zeros(m,n);
      end
      parfor dim = 1:n
         xP             = x;
         % Change Prob to ProbP in rEval call if using the 3 next lines
         %ProbP         = Prob;
         %ProbP.FDVar   = dim;
         %ProbP.cols    = dim;
         z              = x(dim);
         h              = JTol(dim)*(1+abs(z));
         xP(dim)        = z + h;
         rx_ph          = rEval(rxFunc,xP,Prob,NLP,varargin{:});
         J(:,dim)       = (rx_ph-rx)/h;
      end
   else
      % Case 3. Pattern defined
      [ix,iy,iv]        = find(Pattern);
      iv                = double(iv);
      iD                = unique(iy);
      niD               = length(iD);
      if Prob.LargeScale
         JP             = sparse(m,niD);
      else
         JP             = zeros(m,niD);
      end
      parfor i = 1:niD
         xP             = x;
         dim            = iD(i);
         %iz            = find(dim==iy);
         %ixz           = ix(iz);
         % Change Prob to ProbP in rEval call if using the 4 next lines
         %ProbP         = Prob;
         %ProbP.rows    = ixz;
         %ProbP.FDVar   = dim;
         %ProbP.cols    = dim;
         z              = x(dim);
         h              = JTol(dim)*(1+abs(z));
         xP(dim)        = z + h;
         rx_ph          = rEval(rxFunc,xP,Prob,NLP,varargin{:});
         if Prob.LargeScale
            JP(:,i)     = sparse((rx_ph-rx)/h);
         else
            JP(:,i)     = (rx_ph-rx)/h;
         end
      end
      if niD == n
         J    = JP;
      else
         if Prob.LargeScale
            J       = sparse(m,n);
         else
            J       = zeros(m,n);
         end
         J(:,iD) = JP;
      end

      % Waste of time, too slow, to compute every element
      %niv               = length(iv);
      %parfor i = 1:niv
      %   xP             = x;
      %   ProbP          = Prob;
      %   m              = ix(i);
      %   n              = iy(i);
      %   ProbP.rows     = m;
      %   ProbP.cols     = n;
      %   z              = xP(n);
      %   h              = JTol(n)*(1+abs(z));
      %   xP(n)          = z + h;
      %   rx_ph          = rEval(rxFunc,xP,ProbP,NLP,varargin{:});
      %   iv(i)          = (rx_ph(m)-rx(m))/h;
      %end
      %J=sparse(ix,iy,iv,m,n);

   end
end

function rx = rEval(rxFunc, x, Prob, NLP, varargin)
switch NLP
case 2
   rx = nlp_c(x, Prob, varargin{:});
case 1
   rx = nlp_r(x , Prob, varargin{:} );
otherwise
   rx = feval(rxFunc, x, Prob, varargin{:});
end
% MODIFICATION LOG
%
% 110723  hkh  Made parfor version out of FDJac
% 110725  hkh  Save time avoiding sending FDVar info to nlp_c/nlp_r
