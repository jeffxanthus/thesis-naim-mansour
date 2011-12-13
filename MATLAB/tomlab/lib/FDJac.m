% FDJac
%
% Numerical approximation of the Jacobian of the residuals in nonlinear
% least squares problem or the constraint Jacobian for constraint vector c(x).
%
% function [J,rx] = FDJac(x, Prob, rxFunc, rx, varargin)
%
% rx = r(x) or c(x), i.e. the current residual or constraint vector,
% if available
%
% Implementation based on the algorithm FD in Practical Optimization, page 343.
%
% Sparsity pattern in Prob.JacPattern or Prob.ConsPattern is used.
% Preprocessed efficient sparsity input in Prob.ConIx is utilized

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Nov 23, 1998.  Last modified Aug 13, 2009.

function [J,rx] = FDJac(x, Prob, rxFunc, rx, varargin)

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
      [ix,iy,iv]=find(Pattern);
      iv = double(iv);
      mx = max(ConIx);
      for k = 1:mx
         CI         = find(ConIx==k);
         nCI        = length(CI);
         Prob.FDVar = CI;
         Prob.cols  = CI;
         if length(CI) == 1
            Prob.rows = ix(find(CI==iy));
         else
            %Prob.rows = ix(find(any([ones(length(iy),1)*CI==iy*ones(1,mx)]')));
            Prob.rows = ix(find(any([ones(length(iy),1)*CI==iy*ones(1,nCI)]')));
         end
         z          = x(CI);
         h          = JTol(CI).*(1+abs(z));
         x(CI)      = z + h;
         rx_ph      = rEval(rxFunc,x,Prob,NLP,varargin{:});
         x(CI )     = z;
         for j = 1:length(CI)
             iz     = find(CI(j)==iy);
             iv(iz) = (rx_ph(ix(iz))-rx(ix(iz)))/h(j);
         end
      end
      J=sparse(ix,iy,iv,m,n);
   elseif isempty(Pattern)
      if Prob.LargeScale
         J = sparse(m,n);
      else
         J = zeros(m,n);
      end
      for dim = 1:n
         z          = x(dim);
         Prob.FDVar = dim;
         Prob.cols  = dim;
         h          = JTol(dim)*(1+abs(z));
         x(dim)     = z + h;
         rx_ph      = rEval(rxFunc,x,Prob,NLP,varargin{:});
         x(dim)     = z;
         J(:,dim)   = (rx_ph-rx)/h;
      end
   else
      [ix,iy,iv]=find(Pattern);
      iv = double(iv);
      for dim = 1:n
         iz            = find(dim==iy);
         ixz           = ix(iz);
         Prob.rows     = ixz;
         if ~isempty(iz)
            Prob.FDVar = dim;
            Prob.cols  = dim;
            z          = x(dim);
            h          = JTol(dim)*(1+abs(z));
            x(dim)     = z + h;
            rx_ph      = rEval(rxFunc,x,Prob,NLP,varargin{:});
            x(dim)     = z;
            iv(iz)     = (rx_ph(ixz)-rx(ixz))/h;
         end
      end
      J=sparse(ix,iy,iv,m,n);
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
% 981124  mbk  if rx is empty, J is set empty and return.
% 990306  hkh  Safeguard against x slightly outside bounds
% 990626  hkh  Avoid feval
% 990728  hkh  Wrong output name, r, not rx, from rEval
% 990909  hkh  Handle sparsity patterns (NOT in initial estimate part)
% 000911  hkh  Written more efficiently, new step size choice.
% 020416  hkh  n=Prob.N to allow for longer x for minimax, L1 etc.
% 020423  hkh  Safeguard PHI for 0, if no variation in f(x)
% 031006  ango Send global index FDVar telling index of perturbed variable
% 031214  hkh  Unnecessary J = zeros(m,n); moved to if no Pattern,
% 031214  hkh  Use J = sparse(m.n) if Prob.LargeScale = 1
% 040125  hkh  Set double(iv) to avoid treatment at logical array in Matlab 6.5
% 040304  hkh  Make JTol always be a column vector, if expanded from length 1
% 040331  hkh  Return rx as second parameter
% 040406  hkh  If ConIx defined, use a more efficient method with less calls
% 040406  hkh  Send Prob.FDVar telling indicies of perturbed variables
% 040410  hkh  Send Prob.cols and Prob.rows
% 040413  hkh  Check length of x_L and x_U
% 040505  hkh  Incorrect use of mx as length of CI
% 050124  hkh  Set Prob.JacIx as ConIx if LS, move ConIx=Prob.ConIx definition
% 061117  hkh  More efficient code for safeguard against x outside simple bounds
% 090717  med  residual calculation updated
% 090813  med  mlint check