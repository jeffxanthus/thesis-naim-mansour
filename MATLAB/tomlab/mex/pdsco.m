function [x,y,z,inform,PDitns,CGitns,time] = ...
   pdsco( Fname,Aname,b,m,n,options,x0,y0,z0,xsize,zsize, Prob )

% pdsco.m -- Primal-Dual Barrier Method for Separable Convex Objectives
%--------------------------------------------------------------------------
%------
% [x,y,z,inform,PDitns,CGitns,time] = ...
%    pdsco( 'pdsobj', A      ,b,m,n,options,x0,y0,z0,xsize,zsize );
% [x,y,z,inform,PDitns,CGitns,time] = ...
%    pdsco( 'pdsobj','pdsmat',b,m,n,options,x0,y0,z0,xsize,zsize );
%
% pdsco (Primal-Dual barrier method, Separable Convex Objective)
% solves optimization problems of the form
%
%    minimize    phi(x)  +  1/2 norm(gamma*x)^2  +  1/2 norm(r/delta)^2
%      x,r
%    subject to  Ax + r = b,    x > 0,    r unconstrained.
%
% where
% phi(x) is a smooth separable convex function (possibly linear);
% A      is an m x n matrix, or an operator for forming A*x and A'*y;
% b      is a given m-vector;
% gamma  is a primal regularization parameter (typically small but may be 0);
% delta  is a dual regularization parameter (typically small or 1; must be >0);
% With positive gamma and delta, the primal-dual solution (x,y,z) is
% bounded and unique.
%
% EXTERNAL FUNCTIONS:
% options         = pdscoSet;               provided with pdsco.m
% [obj,grad,hess] = pdsobj( x, Prob );      provided by user
%               y = pdsmat( mode,m,n,x );   provided by user if A isn't explicit
%
% OUTPUT ARGUMENTS:
% x         is the primal solution.
% y         is the dual solution associated with Ax + r = b.
% z         is the dual solution associated with x > 0.
%           x > 0, z > 0 is always true.
%           If (x,y,z,r) is optimal, r = delta^2*y and x.*z should be small.
% inform    = 0 if a solution is found;
%           = 1 if too many iterations were required;
%           = 2 if the linesearch failed too often.
% PDitns    is the number of Primal-Dual Barrier iterations required.
% CGitns    is the number of Conjugate-Gradient iterations required
%           if an iterative solver (LSQR) is used.
% time      is the cpu time used.
%
% INPUT ARGUMENTS:
% 'pdsobj'  defines phi(x) and its gradient and diagonal Hessian.
%           grad and hess are both n-vectors.
%           If phi(x) is the linear function c'x, pdsobj should return
%           [obj,grad,hess] = [c'x, c, zeros(n,1)].
% A         is an explicit m x n matrix (preferably sparse!).
% 'pdsmat'  is used if A is not known explicitly.
%           It returns y = A*x (mode=1) or y = A'*x (mode=2).
% b         is a given right-hand side vector.
% m, n      are the dimensions of A (or the A implied by pdsmat).
% options   is a structure that may be set and altered by pdscoSet.
% options.gamma   defines gamma.  Typical value: gamma = 1.0e-4.
% options.delta   defines delta.  Typical value: delta = 1.0e-4 or 1.0.
%                 Values smaller than 1.0e-6 or 1.0e-7 may affect numerical
%                 reliability.  Increasing delta usually improves convergence.
%                 For least-squares applications, delta = 1.0 is appropriate.
% help pdscoSet   describes the other options.
% x0, y0, z0      provide an initial solution.
% xsize, zsize    are estimates of the biggest x and z at the solution.
%                 They are used to scale (x,y,z).  Good estimates
%                 should improve the performance of the barrier method.
%
% Prob       Tomlab Prob structure. Used as second argument in the
%            call to pdsobj
%            Printing in pdsco if Prob.PriLevOpt > 0 (default = 1)
%
%----------------------------------------------------------------------


% AUTHOR:
%    Michael Saunders, Systems Optimization Laboratory (SOL),
%    Stanford University, Stanford, California, USA.
%    saunders@stanford.edu
%
%
% PRIVATE FUNCTIONS:              GLOBAL VARIABLES:
%    pdsxxxdistrib                   pdsAAA (REMOVED IN TOMLAB)
%    pdsxxxlsqr                      pdsDDD1
%    pdsxxxlsqrmat                   pdsDDD2
%    pdsxxxmat
%
%
% NOTES:
% The matrix A is assumed to be reasonably well scaled: norm(A,inf) =~ 1.
% The vector b and objective phi(x) may be of any size, but be sure that
%    xsize and zsize have sensible values.
% The files defining pdsobj and pdsmat must not be called Fname.m or Aname.m !!!!
%
%
% DEVELOPMENT:
% 20 Jun 1997: Original version, derived from pdlp0.m.
% 19 Mar 1998: gamma, delta, wait added as parameters.
% 29 Mar 1998: LSQR's atol decreased dynamically from atol1 to atol2.
%              It looks like atol1 must be 1.0e-6 or less
%              (and atol2 = 1.0e-8 or less).
% 17 Mar 1998: LSproblem =  1 selects original LS problem for dy.
%              This is best if delta is small.
%              LSproblem =  2 converts it to a damped LS problem.
%              This allows lsqr to work with the smaller operator DA'.
%              It should be safe if delta = 1 say.
% 31 Mar 1998: Found that LSdamp > delta often accelerates convergence.
% 31 Mar 1998: Derived "primal" LS problem, which finds dx before dy.
%              LSproblem = 11 selects "primal" LS problem.
%              LSproblem = 12 converts it to a damped LS problem.
%              Implemented pdsobj, pdsmat, pdsxxx to allow for explicit A.
% 09 May 1998: 2x2 KKT may be symmetrized with X^(1/2).
%              LSproblem = 32 solves it via explicit K2 or SYMMLQ.
%              LSproblem = 33 equilibrates first n cols of K2.
%              pdsyyy.m applies operator K for SYMMLQ.
%              Cond(K2) is less than cond(A*D) -- a good sign.
% 09 May 1998: Initialized mu to be the proverbial mu0*(x'z)/n.
% 14 May 1998: Regard delta as a decreasing parameter like mu
%              to implement the conventional penalty function.
%              LSdamp = max( sqrt(mu), delta )  does the trick.
% 23 May 1998: x0,y0,z0 are now optional input parameters.  z is output.
% 03 Jun 1998: delta0 introduced to allow a choice.
%              delta0 = delta fixes delta throughout (this is the default).
%              delta0 = 1    allows delta to decrease with sqrt(mu).
%              Thus, max( sqrt(mu), delta )  <=  "delta"  <=  delta0.
% 13 Jun 1998: mulast = 0.1*opttol introduced.  No need to let mu get smaller.
% 11 Feb 2000: Added diagonal preconditioning for LSQR when A is explicit
%              (and LSproblem = 1, LSmethod = 3). Improved cond(DA') greatly.
%              Seems to halve the LSQR itns as hoped.
% 12 Feb 2000: Added analogous scaling to rows of A for SYMMLQ (LSproblem = 32)
% 14 Apr 2000: Initialize mu from initial [r; t; Xz], not just from gap.
% 15 Apr 2000: Use P = symmmd(ADDA) on first iteration for later chol(ADDA)
%              (if explicit A and LSmethod = 1).
% 16 Apr 2000: mulast = 0.01*opttol seems necessary on Xiaoming's LP problems.
% 28 Sep 2000: options implemented as a structure.
% 21 Nov 2000: Default x0, z0 scaled by norm(b).  (Altered 18 Jan 2001.)
% 24 Nov 2000: Iteration log prints scaled quantities Pinf, Dinf, Cinf.
% 14 Dec 2000: Removed methods that work with 2x2 and 3x3 systems and SYMMLQ.
% 18 Jan 2001: x, y, z, obj, grad scaled to allow for norm(b).
% 21 Jan 2001: Complementarity test is now  Cinf = maxXz <= opttol.
%              Should be ok because maxXz -> final mu = 0.1*opttol.
% 25 Jan 2001: x0, y0, z0, xsize, zsize now required input.
%              x, b are scaled by xsize.
%              y, z are scaled by zsize.
% 29 Jan 2001: Now that problem is scaled, use absolute value mufirst = mu0.
%              For warm starts, user should probably set mu0 = x0min * z0min
%              to get most variables centered.
% 31 Jan 2001: Set atol = max(Pinf, Dinf, Cinf) * 0.1  as in Inexact Newton.
%              atol2 no longer used.
% 07 Feb 2001: Satellite problem seems to need more accurate solves (.01).
% 09 Feb 2001: Test on r3ratio enforces Inexact Newton condition.
% 12 Feb 2001: lsqr is now private function pdsxxxlsqr.
%              Specialized to reduce atol and continue if necessary.
% 14 Feb 2001: stepx and stepz optimized to reduce the merit function f,
%              using Byunggyoo Kim's algebra.
% 13 Apr 2001: PDitns, CGitns now output parameters.
% 25 Apr 2001: Now that starting conditions are better, revert to
%              0.1 rather than the more conservative (and expensive) 0.01
%              in the Inexact Newton test.  (See r3norm below.)
% 26 Apr 2001: Following the central path closely seems safe but slow.
%              Smaller mu0 seems to help.
%              bigcenter = 1e+3 introduced (again) in place of 100.
% 01 Aug 2005: Changed isstr to ischar
%-----------------------------------------------------------------------

%global pdsAAA pdsDDD1 pdsDDD2 pdsDDD3;

global pdsDDD1 pdsDDD2 pdsDDD3;

if nargin < 12
   Prob.PriLevOpt = 1;
   PriLev = 1;
else
   PriLev=DefPar(Prob,'PriLevOpt',1);
end
wait      = options.wait;
% HKH - Avoid keyboard by putting comments on the keyboard commands
if wait > 0
   PriLev = 1;
end

if PriLev > 0
  fprintf('\n   pdsco.m                               Version of 26 Apr 2001')
  fprintf('\n   Primal-dual barrier algorithm to minimize a separable convex')
  fprintf('\n   objective function subject to constraints Ax + r = b, x >= 0\n')
end

operator =  ischar(Aname) | isa(Aname,'function_handle');
explicit = ~operator;

if explicit             % assume Aname is an explicit matrix A.
  nnzA   = nnz(Aname);
  if PriLev > 0
   if  issparse(Aname)
     fprintf('The matrix A is an explicit sparse matrix')
   else
     fprintf('The matrix A is an explicit dense matrix' )
   end
   fprintf('\n\nm        = %8g     n        = %8g      nnz(A)  =%9g', m,n,nnzA)
  end
else
  if ischar(Aname)
     fname   = Aname;
  else  % assume Aname is a function handle
     fstruct = functions(Aname);
     fname   = fstruct.function;
  end
  if PriLev > 0
     disp(['The matrix A is an operator defined by ' fname])
     fprintf('\nm        = %8g     n        = %8g', m,n)
  end
end

%operator =  ischar(Aname);
%explicit = ~operator;

%if explicit               % Aname is an explicit matrix A.
%   [m,n]  =  size(Aname); % In case the user set (m,n) incorrectly.
%   pdsAAA =  Aname;       % Explicit A = global pdsAAA.
%   Aname  = 'pdsxxxmat';  % Use default Ax, A'y routine pdsxxxmat.m.
%   nnzA   =  nnz(pdsAAA);
%   if PriLev > 0
%     if issparse(pdsAAA)
%        fprintf('\nA is an explicit sparse matrix')
%     else
%        fprintf('\nA is an explicit dense matrix' )
%     end
%     fprintf('\nm        = %8g      n       = %8g       nnz(A)  =%9g', m,n,nnzA)
%   end
%else
%   if PriLev > 0
%     fprintf('\nA is an operator')
%     fprintf('\nm        = %8g      n       = %8g', m,n)
%   end
%end

normb  = norm(b ,inf);   normx0 = max(x0);
normy0 = norm(y0,inf);   normz0 = max(z0);
if PriLev > 0
  fprintf('\nmax |b | = %8g      max x0  = %8g', normb , normx0)
  fprintf(                 '       xsize   = %8g', xsize)
  fprintf('\nmax |y0| = %8g      max z0  = %8g', normy0, normz0)
  fprintf(                 '       zsize   = %8g', zsize)
end

%-----------------------------------------------------------------------
% Initialize.
%-----------------------------------------------------------------------
true   = 1;        em   = ones(m,1);
false  = 0;        en   = ones(n,1);
nb     = n + m;    nkkt = nb;
CGitns = 0;
inform = 0;

%-----------------------------------------------------------------------
% Grab input options.
%-----------------------------------------------------------------------
gamma     = options.gamma;
delta     = options.delta;
maxitn    = options.MaxIter;
featol    = options.FeaTol;
opttol    = options.OptTol;
steptol   = options.StepTol;
x0min     = options.x0min;
z0min     = options.z0min;
mu0       = options.mu0;
LSmethod  = options.Method;   % 1=Cholesky    2=QR    3=LSQR (SOL)
  % 4 = Tomlab Tlsqr, special PDCO interface avoiding any callbacks
LSproblem = options.LSproblem;  % See below
itnlim    = options.LSQRMaxIter * min(m,n);
atol1     = options.LSQRatol1;  % Initial  atol
atol2     = options.LSQRatol2;  % Smallest atol (unless atol1 is smaller)
conlim    = options.LSQRconlim;
%wait      = options.wait;

% LSproblem:        %  1 = dy          2 = dy shifted, DLS
                    % 11 = s          12 =  s shifted, DLS    (dx = Ds)
                    % 21 = dx
                    % 31 = 3x3 system, symmetrized by Z^{1/2}
                    % 32 = 2x2 system, symmetrized by X^{1/2}

%-----------------------------------------------------------------------
% Set other parameters.
%-----------------------------------------------------------------------
kminor    = 0;      % 1 stops after each iteration
delta0    = delta;  % Initial value of delta
eta       = 1e-4;   % Linesearch tolerance for "sufficient descent"
maxf      = 10;     % Linesearch backtrack limit (function evaluations)
maxfail   = 1;      % Linesearch failure limit (consecutive iterations)
bigcenter = 1e+3;   % mu is reduced if center < bigcenter.

%if delta >= 0.1,  LSproblem = 2; end

% Parameters for LSQR.
atolmin   = eps;    % Smallest atol if linesearch back-tracks
btol      = 0;      % Should be small (zero is ok)
show      = false;  % Controls lsqr iteration log

if PriLev > 0
  fprintf('\n\nx0min    = %8g      featol  = %8.1e', x0min, featol)
  fprintf(                 '       gamma   = %8.1e', gamma)
  fprintf(  '\nz0min    = %8g      opttol  = %8.1e', z0min, opttol)
  fprintf(                 '       delta   = %8.1e', delta)
  fprintf(  '\nmu0      = %8.1e      steptol = %8g', mu0  , steptol)
  fprintf(                 '       delta0  = %8.1e', delta0)
  fprintf(  '\nbigcenter= %8g'                     , bigcenter)

  fprintf('\n\nLSQR:')
  fprintf('\natol1    = %8.1e      atol2   = %8.1e', atol1 , atol2 )
  fprintf(                 '       btol    = %8.1e', btol )
  fprintf('\nconlim   = %8.1e      itnlim  = %8g'  , conlim, itnlim)
  fprintf(                 '       show    = %8g'  , show )

  %LSmethod  = 3;  %%% Hardwire LSQR
  %LSproblem = 1;  %%% and LS problem defining "dy".

  %if wait
  %    fprintf('\n\nReview parameters... then type "return"\n')
  %    keyboard
  %end

  if eta < 0
      fprintf('\n\nLinesearch disabled by eta < 0')
  end
end

%------------------------------------------------------------------------
% All parameters have now been set.
%------------------------------------------------------------------------
time    = cputime;
useChol = LSmethod == 1;
useQR   = LSmethod == 2;
direct  = LSmethod <= 2 & explicit;
solver  = '  LSQR';   if LSproblem >= 30, solver  = 'SYMMLQ'; end

%------------------------------------------------------------------------
% Scale the input data.
% The scaled variables are
%    xbar     = x/beta,
%    ybar     = y/zeta,
%    zbar     = z/zeta.
% Define
%    theta    = beta*zeta;
% The scaled function is
%    phibar   = ( 1   /theta) fbar(beta xbar),
%    gradient = (beta /theta) grad,
%    Hessian  = (beta2/theta) hess.
%------------------------------------------------------------------------
beta   = xsize;   if beta==0, beta = 1; end    % beta scales b, x.
zeta   = zsize;   if zeta==0, zeta = 1; end    % zeta scales y, z.
theta  = beta*zeta;                            % theta scales obj.
gamma  = gamma*beta/sqrt(theta);
delta  = delta*zeta/sqrt(theta);

beta2  = beta^2;    gamma2 = gamma^2;
delta2 = delta^2;   delta0 = delta;

b      = b /beta;   y0 = y0/zeta;
x0     = x0/beta;   z0 = z0/zeta;

%------------------------------------------------------------------------
% Initialize x, y, z, objective, etc.
%------------------------------------------------------------------------
x      = max( x0, x0min );   clear x0
z      = max( z0, z0min );   clear z0
y      = y0;                 clear y0

[obj,grad,hess] = feval( Fname, (x*beta), Prob );
obj    = obj        /theta;               % Scaled obj.
grad   = grad*(beta /theta) + gamma2*x;   % grad includes x regularization.
Q      = hess*(beta2/theta) + gamma2;     % Q    includes x regularization.

%------------------------------------------------------------------------
% Initialize linear residuals  rlin =   b - A*x,
%                              tlin = - z - A'*y.
%------------------------------------------------------------------------
if explicit
   rlin =   b - Aname *x;
   tlin = - z - Aname'*y;
else
   rlin =   b - pdsxxxmat( Aname, 1, m, n, x );
   tlin = - z - pdsxxxmat( Aname, 2, m, n, y );
end

%------------------------------------------------------------------------
% Initialize mu
% and the nonlinear residuals  r    = rlin - delta^2*y,
%                              t    = tlin + gamma^2*x + grad,
%                              v    = mu*e - X*z.
%
% 09 May 1998: Initialize mu to a fraction of x'z/n (as others do).
% 14 May 1998: Make delta decrease with mu.
% 03 Jun 1998: Keep delta <= delta0.
% 13 Jun 1998: Keep mu >= mulast.
% 14 Apr 2000: mufirst = mu0 * (x'z)/n is too small if initial r, t are big.
%              Revert to balancing two parts of rhs of KKT system:
%              (  b - Ax   )     ( 0  )     (  r  )
%              (g - A'y - z)  =  ( 0  )  +  (  t  ).
%              (mu e  -  Xz)     (mu e)     (- Xz )
%
% 25 Jan 2001: Now that b and obj are scaled (and hence x,y,z),
%              we should be able to use an absolute value: mufirst = mu0;
%              But 0.1 worked poorly on StarTest1 (with x0min = z0min = 0.1).
%
% 29 Jan 2001: We might as well use mu0 = x0min * z0min;
%              so that most variables are centered after a warm start.
%------------------------------------------------------------------------
r       = rlin - delta2*y;
t       = tlin + grad;           % grad includes gamma2*x
Xz      = x.*z;

%f      = norm([r; t; Xz]);      % Original norm.
%f      = norm([r; t; Xz],1);    % Probably should have been that.
%mufirst= mu0 * f / n;           % Replaces  mufirst = mu0 * (sum(Xz) / n);

mufirst = mu0;                   % Absolute value.
mulast  = 0.1 * opttol;
mulast  = min( mulast, mufirst );
mu      = mufirst;
LSdamp  = max( sqrt(mu), delta );
LSdamp  = min( LSdamp  , delta0);
LSdamp2 = LSdamp^2;
r       = rlin - LSdamp2*y;      % 25 Jan 2001: Should have been there earlier.
v       = mu   - Xz;
f       = norm([r; t; v]);       % f = merit function for the linesearch.

% Initialize other things.

PDitns    = 0;
converged = 0;
atol      = atol1;
atol2     = max( atol2, atolmin );

%  Iteration log.

stepx   = 0;
stepz   = 0;
nf      = 0;
itncg   = 0;
nfail   = 0;

if PriLev > 0
   head1   = '\n\nItn   mu stepx stepz  Pinf  Dinf';
   head2   = '  Cinf   Objective    nf  center';
   if direct, head3 = ''; else head3 =['  atol ' solver ' Inexact']; end
   fprintf([ head1 head2 head3 ])
end

mininf  = 1e-99;
regterm = gamma2 * (x'*x)  +  LSdamp2 * (y'*y);
objreg  = obj  +  0.5 * regterm;
objtrue = objreg * theta;

maxXz   = max(Xz);
minXz   = min(Xz);
center  = maxXz / minXz;

Pinf    = norm(r,inf);   Pinf = max( Pinf, mininf );
Dinf    = norm(t,inf);   Dinf = max( Dinf, mininf );
Cinf    = maxXz;         Cinf = max( Cinf, mininf );

if PriLev > 0
   fprintf('\n%3g                 ', PDitns       )
   fprintf('%6.1f%6.1f' , log10(Pinf), log10(Dinf))
   fprintf('%6.1f%15.7e', log10(Cinf), objtrue    )
   fprintf('   %8.1f'   , center                  )

   if kminor
      fprintf('\n\nStart of first minor itn...\n')
      %keyboard
   end
end

%-----------------------------------------------------------------------
% Main loop.
%-----------------------------------------------------------------------

while ~converged
   PDitns = PDitns + 1;
%  x1 = x;  y1 = y;  z1 = z;  % Save for debugging at end of iteration.
%  t1 = t;  r1 = r;  v1 = v;

% 31 Jan 2001: Set atol according to actual progress, a la Inexact Newton.
% 07 Feb 2001: 0.1 not small enough for Satellite problem.  Try 0.01.
% 25 Apr 2001: 0.01 seems wasteful for Star problem.
%              Now that starting conditions are better, go back to 0.1.

   r3norm = max([Pinf  Dinf  Cinf]);
   atol   = min([atol  r3norm*0.1]);
   atol   = max([atol  atolmin   ]);

%-----------------------------------------------------------------------
%  Define a damped Newton iteration for solving
%     ( r, t, v ) = 0,
%  keeping  x, z > 0.  We eliminate dz
%  to obtain the system
%
%     [-H   A'] [ dx ] = [ w ],    d2I = delta^2 I,   w = t - v./x;
%     [ A  d2I] [ dy ] = [ r ]
%
%  which is equivalent to
%
%     [-dI  DA'] [ s  ] = [  Dw   ],    dI = delta I,
%     [ AD  dI ] [ dy ] = [r/delta]     D  = H^{-1/2),  dx = delta D s,
%
%  and also equivalent to the least-squares problem
%
%     min | [ DA']dy  -  [  Dw   ] |.                                (*)
%         | [ dI ]       [r/delta] |
%
% 17 Mar 1998: We solve the latter as the damped least-squares problem
%
%     min | [ DA']dybar  -  [D wbar] |,  wbar = w - (A'r)/delta^2,   (**)
%         | [ dI ]          [  0   ] |     dy = dybar + r/delta^2,
%
%    to allow lsqr to work with the smaller operator DA'.
%
% 30 Mar 1998: LSproblem = 1 or 2 selects (*) or (**) respectively.
%              (*)  seems safer if delta is small, but
%              (**) is more efficient (and safe) if delta = 1, say.
%
% 31 Mar 1998: LSproblem = 11 or 12 selects alternative LS problem
%
%     min | [ AD ]s  -  [    r    ] |,     dx = Ds,                 (***)
%         | [ dI ]      [-delta Dw] |      dy = (r - Adx) / delta^2
%
% and associated damped LS problem:
%
%     min | [ AD ]sbar  -  [ rbar ] |,      s = sbar - Dw,         (****)
%         | [ dI ]         [  0   ] |    rbar = r + AD^2 w.
%
% 06 Apr 1998: LSproblem = 21 selects equivalent LS problem
%
%     min | [    A     ]dx  -  [    r    ] |,                     (*****)
%         | [delta Dinv]       [-delta Dw] |     dy = (r - Adx) / delta^2
%-----------------------------------------------------------------------
   H       = Q  +  z./x;
   w       = t  -  v./x;
   Hinv    = 1 ./ H;
   D       = sqrt(Hinv);
   pdsDDD1 = D;
   rw      = [explicit  LSproblem  LSmethod  LSdamp m n 0];  % Passed to LSQR.

   if LSproblem == 1
      %-----------------------------------------------------------------
      %  Solve (*) for dy.
      %-----------------------------------------------------------------
      Dw   = D.*w;

      if direct
         DD   = sparse( 1:n, 1:n, D, n, n );
         AD   = Aname*DD;

         if useChol
            d2I  = LSdamp2 * em;
            d2I  = sparse( 1:m, 1:m, d2I, m, m );
            ADDA = AD * AD'  +  d2I;

            if PDitns==1,  P = symmmd(ADDA);  end  % Do ordering only once.

            [R,indef] = chol(ADDA(P,P));
            if indef
               fprintf('\n\n   chol says AD^2A'' is not positive definite')
               break
            end

            rhs  = r  +  Aname * (Hinv.*w);
%           dy   = ADDA \ rhs;
            dy   = R \ (R' \ rhs(P));   dy(P) = dy;
         
	 else % useQR
            dI   = LSdamp * em;
            rhs  = [ Dw; r/LSdamp ];
            dy   = [ AD'; diag(dI) ] \ rhs;
         end

         Ady  = Aname'*dy;
         dx   = Hinv .* (Ady - w);
         Adx  = Aname*dx;


      else % Iterative solve using LSQR.
         rhs     = [ Dw; r/LSdamp ];
         damp    = 0;

         if explicit            % A is a sparse matrix.
	    precon     = true;
	    if precon    % Construct diagonal preconditioner for LSQR
	       DD      = sparse( 1:n, 1:n, D, n, n );
	       AD      = Aname*DD;
	       AD2     = AD.^2;
	       wD      = sum( (AD2') )';   % Sum of squares of each row of AD
	       wD      = sqrt( wD + LSdamp2 );
	       pdsDDD2 = 1 ./ wD;
	    end
         else                   % A is an operator
	    precon    = false;  
	 end

	 rw(7)        = precon;
         info.atolmin = atolmin;
	 info.r3norm  = f;    % Must be the 2-norm here.

% HKH - Add options to call Tomlab
        switch LSmethod
        case 3
	  [ dy, istop, itncg, outfo ] = ...
	      pdsxxxlsqr( nb, m, 'pdsxxxlsqrmat', Aname, rw, rhs, damp, ...
 	                 atol, btol, conlim, itnlim, show, info );

	  if precon, dy = pdsDDD2 .* dy; end
        otherwise % Default use special PDCO LSQR MEX interface in Tomlab
          [ dy, istop, itncg, outfo ] = ...
              pdsTlsqr( nb,m,'pdsTlsqrmat',Aname,rw,rhs,damp, ...
                       atol,btol,conlim,itnlim,show,info );
          % Done in special version of LSQR: if precon, dy = pdsDDD2 .* dy; end
        end


	 if istop == 3 | istop == 7   % conlim or itnlim
	    fprintf('\n    LSQR stopped early:  istop = %3d', istop)
	 end

	 atolold = outfo.atolold;
	 atol    = outfo.atolnew;
	 r3ratio = outfo.r3ratio;

         Ady = pdsxxxmat( Aname, 2, m, n, dy );   % A'dy
	 dx  = Hinv .* (Ady - w);             %   dx
	 Adx = pdsxxxmat( Aname, 1, m, n, dx );   % A dx
      end % LSproblem 1
   
   else
      disp( 'This LSproblem not yet implemented' )
      disp( 'Major failure in pdsco - contact Tomlab Optimization:' )
      disp( '   support@tomopt.com' )
      error('pdsco: Illegal LSproblem!!!')
      %keyboard
   end
%-----------------------------------------------------------------------

   CGitns = CGitns + itncg;

%-----------------------------------------------------------------------
%  dx and dy are now known.  Get dz.
%  dz  = xinv.*(v - z.*dx);  % is the classical formula -- no good???
%-----------------------------------------------------------------------
   if LSproblem ~= 31
      dz  = t - Ady + Q.*dx;
   end

%-----------------------------------------------------------------------
%   Find the maximum step.
%-----------------------------------------------------------------------
   stepx = 1.0e+20;
   stepz = 1.0e+20;

   blocking  =  find( dx < 0 );
   if length( blocking ) > 0
      steps  =  x(blocking) ./ (- dx(blocking));
      stepx  =  min( steps );
   end

   blocking  =  find( dz < 0 );
   if length( blocking ) > 0
      steps  =  z(blocking) ./ (- dz(blocking));
      stepz  =  min( steps );
   end

   stepxmax =  stepx;
   stepzmax =  stepz;
   stepx    =  min( steptol * stepx, 1 );
   stepz    =  min( steptol * stepz, 1 );

%-----------------------------------------------------------------------
%  Optimize stepx and stepz (Byung's 1-D search).
%-----------------------------------------------------------------------
   optsteps = 0;
   if optsteps
    if (stepx > stepz)
      normAdx2   = Adx'*Adx;
      rtAdx      = r'*Adx;
      dxz        = dx .* z;
      normdxz2   = dxz'*dxz;
      gammadx    = gamma2 * dx;
      gammadx2   = gammadx' * gammadx;
      vtdxz      = v'*dxz;
      ttgammadx  = t'*gammadx;
      numerator  = rtAdx + vtdxz - ttgammadx;
      denominat  = normAdx2 + normdxz2 + gammadx2;
      alphax     = (1 - stepz) * numerator / denominat;
      stepxmax   = stepz + alphax;
      stepxmax   = max(stepxmax, stepz);
      stepx      = min(stepxmax, stepx);
    elseif (stepx < stepz)
      normdy2    = dy'*dy;
      rtdy       = r'*dy;
      xdz        = x .* dz;
      normxdz2   = xdz'*xdz;
      vtxdz      = v'*xdz;
      Adydz      = Ady + dz;
      ttAdydz    = t'*Adydz;
      Adydz2     = Adydz'*Adydz;
      numerator  = rtdy + vtxdz + ttAdydz;
      denominat  = normdy2 + normxdz2 + Adydz2;
      alphaz     = (1 - stepx) * numerator / denominat;
      stepzmax   = stepx + alphaz;
      stepzmax   = max(stepzmax, stepx);
      stepz      = min(stepzmax, stepz);
    end
   end % Byung's steps

%-----------------------------------------------------------------------
%  Backtracking linesearch.
%-----------------------------------------------------------------------
   fail     =  true;
   nf       =  0;

   while nf < maxf
      nf       =   nf + 1;
      xnew     =   x  +  stepx * dx;
      ynew     =   y  +  stepz * dy;
      znew     =   z  +  stepz * dz;
      [obj,grad,hess] = feval( Fname, (xnew*beta), Prob );
      obj      =   obj /theta;
      grad     =   grad*(beta /theta)  +  gamma2*xnew;
      Q        =   hess*(beta2/theta)  +  gamma2;
      rlinnew  =   rlin     -  stepx * Adx;
      tlinnew  =   tlin     -  stepz *(Ady + dz);
      rnew     =   rlinnew  -  LSdamp2 * ynew;
      tnew     =   tlinnew  +  grad;
      Xznew    =   xnew .* znew;
      vnew     =   mu  -  Xznew;
      fnew     =   norm([rnew;  tnew;  vnew]);     % Must be 2-norm here.
      step     =   min( stepx, stepz );
      % if fnew >= f,
      %    fratio = fnew / f;
      %    fprintf( 'fnew / f%8.1e', fratio)
      % end
      if fnew <= (1 - eta*step)*f
         fail = false;
         break;
      end

%     If the first attempt doesn't work,
%     make stepx and stepz the same.

      if nf == 1 & stepx ~= stepz
         stepx   = min( stepx, stepz );
         stepz   = stepx;
      elseif nf < maxf
         stepx   = stepx/2;
         stepz   = stepx;
      end;
   end

   if fail
      fprintf('\n     Linesearch failed (nf too big)');
      nfail = nfail + 1;
   else
      nfail = 0;
   end

   x = xnew;   r = rnew;   rlin = rlinnew;   f  = fnew;
   y = ynew;   t = tnew;   tlin = tlinnew;   Xz = Xznew;
   z = znew;   v = vnew;

%-----------------------------------------------------------------------
%  Set convergence measures.
%-----------------------------------------------------------------------

   regterm = gamma2 * (x'*x)  +  LSdamp2 * (y'*y);
   objreg  = obj  +  0.5*regterm;
   objtrue = objreg * theta;

   maxXz   = max(Xz);
   minXz   = min(Xz);
   center  = maxXz / minXz;

   Pinf    = norm(r,inf);   Pinf = max( Pinf, mininf );
   Dinf    = norm(t,inf);   Dinf = max( Dinf, mininf );
   Cinf    = maxXz;         Cinf = max( Cinf, mininf );
   primalfeas    = Pinf  <=  featol;
   dualfeas      = Dinf  <=  featol;
   complementary = Cinf  <=  opttol;
   enough        = PDitns>=       4;  % Prevent premature termination.
   converged     = primalfeas  &  dualfeas  &  complementary  &  enough;

%-----------------------------------------------------------------------
%  Iteration log.
%-----------------------------------------------------------------------
if PriLev > 0
   str1    = sprintf('\n%3g%5.1f' , PDitns     , log10(mu)   );
   str2    = sprintf('%6.3f%6.3f' , stepx      , stepz       );
   if stepx < 0.001 | stepz < 0.001
      str2 = sprintf('%6.1e%6.1e' , stepx      , stepz       );
   end

   str3    = sprintf('%6.1f%6.1f' , log10(Pinf), log10(Dinf));
   str4    = sprintf('%6.1f%15.7e', log10(Cinf), objtrue     );
   str5    = sprintf('%3g%8.1f'   , nf         , center      );
   if center > 99999
      str5 = sprintf('%3g%8.1e'   , nf         , center      );
   end
   fprintf([str1 str2 str3 str4 str5])
   if direct
      % relax
   else
      fprintf(' %5.1f%7g%7.3f', log10(atolold), itncg, r3ratio)
   end

%-----------------------------------------------------------------------
%  Test for termination.
%-----------------------------------------------------------------------
   if kminor
      fprintf( 'Start of next minor itn...\n')
      %keyboard
   end
end

   if converged
      if PriLev > 0
         fprintf('\n   Converged')
      end
   elseif PDitns >= maxitn
      if PriLev > 0
         fprintf('\n   Too many iterations')
      end
      inform = 1;
      break
   elseif nfail  >= maxfail
      if PriLev > 0
         fprintf('\n   Too many linesearch failures')
      end
      inform = 2;
      break
   else

      % Reduce mu and LSdamp, and reset certain residuals.

      stepmu  = min( stepx , stepz   );
      stepmu  = min( stepmu, steptol );
      muold   = mu;
      mu      = mu   -  stepmu * mu;
      if center >= bigcenter,  mu = muold;  end

%     mutrad  = mu0 * (sum(Xz) / n); % 24 May 1998:  Traditional value. But--
%     mu      = min( mu, mutrad );   % it seemed to decrease mu too much.

      mu      = max( mu, mulast );   % 13 Jun 1998:  No need for smaller mu.
      LSdamp  = max( sqrt(mu), delta );
      LSdamp  = min( LSdamp  , delta0);
      LSdamp2 = LSdamp^2;
      r       = rlin  -  LSdamp2 * y;
      v       = mu  -  Xz;
      f       = norm([r; t; v]);

%     Reduce atol for LSQR (and SYMMLQ).
%     NOW DONE AT TOP OF LOOP.

      atolold = atol;
      % if atol > atol2
      %   atolfac = (mu/mufirst)^0.25;
      %   atol    = max( atol*atolfac, atol2 );
      % end

%     atol    = min( atol, mu );     % 22 Jan 2001:  a la Inexact Newton.
%     atol    = min( atol, 0.5*mu ); % 30 Jan 2001:  A bit tighter

      % If the linesearch took more than one function (nf > 1),
      % we assume the search direction needed more accuracy
      % (though this may be true only for LPs).
      % 12 Jun 1998: Ask for more accuracy if nf > 2.
      % 24 Nov 2000: Also if the steps are small.
      % 30 Jan 2001: Small steps might be ok with warm start.
      % 06 Feb 2001: Not necessarily.  Reinstated tests in next line.

      if nf > 2 | min( stepx, stepz ) <= 0.01
         atol = atolold*0.1;
      end
   end
end
%-----------------------------------------------------------------------
% End of main loop.
%-----------------------------------------------------------------------

if PriLev > 0
   fprintf('\n\nmax( x ) =%10.3f', max(x)     )
   fprintf('    max(|y|) =%10.3f', max(abs(y)))
   fprintf('    max( z ) =%10.3f', max(z)     )
   fprintf('   scaled')
end

x = x*beta;   y = y*zeta;   z = z*zeta;   % Unscale x, y, z.

if PriLev > 0
   fprintf(  '\nmax( x ) =%10.3f', max(x)     )
   fprintf('    max(|y|) =%10.3f', max(abs(y)))
   fprintf('    max( z ) =%10.3f', max(z)     )
   fprintf(' unscaled')
end

time   = cputime - time;

if PriLev > 0
   str1   = sprintf('\nPDitns   =%10g', PDitns );
   str2   = sprintf(     ' itns =%10g', CGitns );
   fprintf( [str1 ' ' solver str2] )
   fprintf('    time     =%10.1f', time);

   pdsxxxdistrib(x,z);   % Private function

   %if wait
   %   keyboard
   %end
end
%-----------------------------------------------------------------------
% End function pdsco.m
%-----------------------------------------------------------------------


function pdsxxxdistrib(x,z)

% pdsxxxdistrib(x) or pdsxxxdistrib(x,z) prints the
% distribution of 1 or 2 vectors.
%
% 18 Dec 2000.  First version with 2 vectors.

two  = nargin > 1;
fprintf('\n\nDistribution of vector     x')
if two, fprintf('         z'); end

x1   = 10^(floor(log10(max(x))) + 1);
z1   = 10^(floor(log10(max(z))) + 1);
x1   = max(x1,z1);
kmax = 10;

for k = 1:kmax
    x2 = x1;    x1 = x1/10;
    if k==kmax, x1 = 0; end
    nx = length(find(x>=x1 & x<x2));
    fprintf('\n[%7.3g,%7.3g )%10g', x1, x2, nx);
    if two
       nz = length(find(z>=x1 & z<x2));
       fprintf('%10g', nz);
    end
end

disp(' ')

%-----------------------------------------------------------------------
% End private function pdsxxxdistrib
%-----------------------------------------------------------------------


function [ x, istop, itn, outfo ] = ...
   pdsxxxlsqr( m, n, aprodname, iw, rw, b, damp, ...
	       atol, btol, conlim, itnlim, show, info )

% Special version of LSQR for use with pdsco.m.
% It continues with a reduced atol if a pdsco-specific test isn't
% satisfied with the input atol.
%
% LSQR solves  Ax = b  or  min |b - Ax|_2  if damp = 0,
% or   min | (b)  -  (  A   )x |   otherwise.
%          | (0)     (damp I)  |2
% A  is an m by n matrix defined by  y = aprod( mode,m,n,x,iw,rw ),
% where the parameter 'aprodname' refers to a function 'aprod' that
% performs the matrix-vector operations.
% If mode = 1,   aprod  must return  y = Ax   without altering x.
% If mode = 2,   aprod  must return  y = A'x  without altering x.
% WARNING:   The file containing the function 'aprod'
%            must not be called aprodname.m !!!!

%-----------------------------------------------------------------------
% LSQR uses an iterative (conjugate-gradient-like) method.
% For further information, see 
% 1. C. C. Paige and M. A. Saunders (1982a).
%    LSQR: An algorithm for sparse linear equations and sparse least squares,
%    ACM TOMS 8(1), 43-71.
% 2. C. C. Paige and M. A. Saunders (1982b).
%    Algorithm 583.  LSQR: Sparse linear equations and least squares problems,
%    ACM TOMS 8(2), 195-209.
% 3. M. A. Saunders (1995).  Solution of sparse rectangular systems using
%    LSQR and CRAIG, BIT 35, 588-604.
%
% Input parameters:
% iw, rw      are not used by lsqr, but are passed to aprod.
% atol, btol  are stopping tolerances.  If both are 1.0e-9 (say),
%             the final residual norm should be accurate to about 9 digits.
%             (The final x will usually have fewer correct digits,
%             depending on cond(A) and the size of damp.)
% conlim      is also a stopping tolerance.  lsqr terminates if an estimate
%             of cond(A) exceeds conlim.  For compatible systems Ax = b,
%             conlim could be as large as 1.0e+12 (say).  For least-squares
%             problems, conlim should be less than 1.0e+8.
%             Maximum precision can be obtained by setting
%             atol = btol = conlim = zero, but the number of iterations
%             may then be excessive.
% itnlim      is an explicit limit on iterations (for safety).
% show = 1    gives an iteration log,
% show = 0    suppresses output.
% info        is a structure special to pdsco.m, used to test if
%             was small enough, and continuing if necessary with smaller atol.
%
%
% Output parameters:
% x           is the final solution.
% istop       gives the reason for termination.
% istop       = 1 means x is an approximate solution to Ax = b.
%             = 2 means x approximately solves the least-squares problem.
% rnorm       = norm(r) if damp = 0, where r = b - Ax,
%             = sqrt( norm(r)**2  +  damp**2 * norm(x)**2 ) otherwise.
% xnorm       = norm(x).
% var         estimates diag( inv(A'A) ).  Omitted in this special version.
% outfo       is a structure special to pdsco.m, returning information
%             about whether atol had to be reduced.
%             
% Other potential output parameters:
% anorm, acond, arnorm, xnorm
%
%        1990: Derived from Fortran 77 version of LSQR.
% 22 May 1992: bbnorm was used incorrectly.  Replaced by anorm.
% 26 Oct 1992: More input and output parameters added.
% 01 Sep 1994: Matrix-vector routine is now a parameter 'aprodname'.
%              Print log reformatted.
% 14 Jun 1997: show  added to allow printing or not.
% 30 Jun 1997: var   added as an optional output parameter.
%              It returns an estimate of diag( inv(A'A) ).
% 12 Feb 2001: atol  can now be reduced and iterations continued if necessary.
%              info, outfo are new problem-dependent parameters for such purposes.
%              In this version they are specialized for pdsco.m.
%-----------------------------------------------------------------------

%     Initialize.

msg=['The exact solution is  x = 0                              '
     'Ax - b is small enough, given atol, btol                  '
     'The least-squares solution is good enough, given atol     '
     'The estimate of cond(Abar) has exceeded conlim            '
     'Ax - b is small enough for this machine                   '
     'The least-squares solution is good enough for this machine'
     'Cond(Abar) seems to be too large for this machine         '
     'The iteration limit has been reached                      '];

%wantvar= nargout >= 6;
%if wantvar, var = zeros(n,1); end
  
itn    = 0;		istop  = 0;		nstop  = 0;
ctol   = 0;		if conlim > 0, ctol = 1/conlim; end;
anorm  = 0;		acond  = 0;
dampsq = damp^2;	ddnorm = 0;		res2   = 0;
xnorm  = 0;		xxnorm = 0;		z      = 0;
cs2    = -1;		sn2    = 0;

% Set up the first vectors u and v for the bidiagonalization.
% These satisfy  beta*u = b,  alfa*v = A'u.

u      = b(1:m);	x    = zeros(n,1);
alfa   = 0;		beta = norm( u );
if beta > 0
   u = (1/beta) * u;	v = feval( aprodname, 2, m, n, u, iw, rw );
   alfa = norm( v );
end
if alfa > 0
   v = (1/alfa) * v;    w = v;
end

arnorm = alfa * beta;   if arnorm == 0, disp(msg(1,:)); return, end

rhobar = alfa;		phibar = beta;		bnorm  = beta;
rnorm  = beta;
head1  = '   Itn      x(1)      Function';
head2  = ' Compatible   LS      Norm A   Cond A';

if show
   disp(' ')
   disp([head1 head2])
   test1  = 1;		test2  = alfa / beta;
   str1   = sprintf( '%6g %12.5e %10.3e',   itn, x(1), rnorm );
   str2   = sprintf( '  %8.1e %8.1e',       test1, test2 );
   disp([str1 str2])
end

%------------------------------------------------------------------
%     Main iteration loop.
%------------------------------------------------------------------
while itn < itnlim
      itn = itn + 1;
%     Perform the next step of the bidiagonalization to obtain the
%     next  beta, u, alfa, v.  These satisfy the relations
%                beta*u  =  a*v   -  alfa*u,
%                alfa*v  =  A'*u  -  beta*v.

      u    = feval( aprodname, 1, m, n, v, iw, rw )  -  alfa*u;
      beta = norm( u );
      if beta > 0
         u     = (1/beta) * u;
         anorm = norm([anorm alfa beta damp]);
         v     = feval( aprodname, 2, m, n, u, iw, rw )  -  beta*v;
         alfa  = norm( v );
         if alfa > 0,  v = (1/alfa) * v; end
      end

%     Use a plane rotation to eliminate the damping parameter.
%     This alters the diagonal (rhobar) of the lower-bidiagonal matrix.

      rhobar1 = norm([rhobar damp]);
      cs1     = rhobar / rhobar1;
      sn1     = damp   / rhobar1;
      psi     = sn1 * phibar;
      phibar  = cs1 * phibar;

%     Use a plane rotation to eliminate the subdiagonal element (beta)
%     of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.

      rho     = norm([rhobar1 beta]);
      cs      =   rhobar1/ rho;
      sn      =   beta   / rho;
      theta   =   sn * alfa;
      rhobar  = - cs * alfa;
      phi     =   cs * phibar;
      phibar  =   sn * phibar;
      tau     =   sn * phi;

%     Update x and w.

      t1      =   phi  /rho;
      t2      = - theta/rho;
      dk      =   (1/rho)*w;

      x       = x      +  t1*w;
      w       = v      +  t2*w;
      ddnorm  = ddnorm +  norm(dk)^2;
%     if wantvar, var = var  +  dk.*dk; end

%     Use a plane rotation on the right to eliminate the
%     super-diagonal element (theta) of the upper-bidiagonal matrix.
%     Then use the result to estimate  norm(x).

      delta   =   sn2 * rho;
      gambar  = - cs2 * rho;
      rhs     =   phi  -  delta * z;
      zbar    =   rhs / gambar;
      xnorm   =   sqrt(xxnorm + zbar^2);
      gamma   =   norm([gambar theta]);
      cs2     =   gambar / gamma;
      sn2     =   theta  / gamma;
      z       =   rhs    / gamma;
      xxnorm  =   xxnorm  +  z^2;

%     Test for convergence.
%     First, estimate the condition of the matrix  Abar,
%     and the norms of  rbar  and  Abar'rbar.

      acond   =   anorm * sqrt( ddnorm );
      res1    =   phibar^2;
      res2    =   res2  +  psi^2;
      rnorm   =   sqrt( res1 + res2 );
      arnorm  =   alfa * abs( tau );

%     Now use these norms to estimate certain other quantities,
%     some of which will be small near a solution.

      test1   =   rnorm / bnorm;
      test2   =   arnorm/( anorm * rnorm );
      test3   =       1 / acond;
      t1      =   test1 / (1    +  anorm * xnorm / bnorm);
      rtol    =   btol  +  atol *  anorm * xnorm / bnorm;

%     The following tests guard against extremely small values of
%     atol, btol  or  ctol.  (The user may have set any or all of
%     the parameters  atol, btol, conlim  to 0.)
%     The effect is equivalent to the normal tests using
%     atol = eps,  btol = eps,  conlim = 1/eps.

      if itn >= itnlim,   istop = 7; end
      if 1 + test3  <= 1, istop = 6; end
      if 1 + test2  <= 1, istop = 5; end
      if 1 + t1     <= 1, istop = 4; end

%     Allow for tolerances set by the user.

      if  test3 <= ctol,  istop = 3; end
      if  test2 <= atol,  istop = 2; end
      if  test1 <= rtol,  istop = 1; end

%-----------------------------------------------------------------------
%     SPECIAL TEST THAT DEPENDS ON pdsco.m.
%     Aname in pdsco  is  iw in lsqr.
%     dy              is  x
%     Other stuff     is in info.
%     We allow for diagonal preconditioning in pdsDDD2.
%-----------------------------------------------------------------------
      if istop > 0
	 r3new     = arnorm;
	 r3ratio   = r3new / info.r3norm;
         atolold   = atol;
	 atolnew   = atol;
	 
	 if atol > info.atolmin
	    if     r3ratio <= 0.1     % dy seems good
	       % Relax
	    elseif r3ratio <= 0.5     % Accept dy but make next one more accurate.
	       atolnew = atolnew * 0.1;
	    else                      % Recompute dy more accurately
	       fprintf('\n                                ')
	       fprintf('                                ')
	       fprintf(' %5.1f%7g%7.3f', log10(atolold), itn, r3ratio)
	       atol    = atol * 0.1;
	       atolnew = atol;
	       istop   = 0;
	    end
	 end

	 outfo.atolold = atolold;
	 outfo.atolnew = atolnew;
      	 outfo.r3ratio = r3ratio;
      end
      
%-----------------------------------------------------------------------
%     See if it is time to print something.
%-----------------------------------------------------------------------
      prnt = 0;
      if n     <= 40       , prnt = 1; end
      if itn   <= 10       , prnt = 1; end
      if itn   >= itnlim-10, prnt = 1; end
      if rem(itn,10) == 0  , prnt = 1; end
      if test3 <=  2*ctol  , prnt = 1; end
      if test2 <= 10*atol  , prnt = 1; end
      if test1 <= 10*rtol  , prnt = 1; end
      if istop ~=  0       , prnt = 1; end

      if prnt == 1
         if show
            str1 = sprintf( '%6g %12.5e %10.3e',   itn, x(1), rnorm );
            str2 = sprintf( '  %8.1e %8.1e',       test1, test2 );
            str3 = sprintf( ' %8.1e %8.1e',        anorm, acond );
            disp([str1 str2 str3])
         end
      end
      if istop > 0, break, end
end

%     End of iteration loop.
%     Print the stopping condition.

if show
   disp(' ')
   disp('LSQR finished')
   disp(msg(istop+1,:))
   disp(' ')
   str1 = sprintf( 'istop  =%8g   itn    =%8g',      istop, itn    );
   str2 = sprintf( 'anorm  =%8.1e   acond  =%8.1e',  anorm, acond  );
   str3 = sprintf( 'rnorm  =%8.1e   arnorm =%8.1e',  rnorm, arnorm );
   str4 = sprintf( 'bnorm  =%8.1e   xnorm  =%8.1e',  bnorm, xnorm  );
   disp([str1 '   ' str2])
   disp([str3 '   ' str4])
   disp(' ')
end

%-----------------------------------------------------------------------
% End private function pdsxxxlsqr
%-----------------------------------------------------------------------


function y = pdsxxxlsqrmat( mode, mlsqr, nlsqr, x, Aname, rw )
%
% pdsxxxlsqrmat is required by pdsco.m (when it calls pdsxxxlsqr.m).
% It forms Mx or M'x for some operator M that depends on LSproblem below.
%
% mlsqr, nlsqr  are the dimensions of the LS problem that lsqr is solving.
%
% Aname is the name of the user's (Ax,A'y) routine
% or the default routine 'pdsxxxmat'.
%
% rw contains parameters [explicit LSproblem LSmethod LSdamp]
% from pdsco.m to say which least-squares subproblem is being solved.
%
% global pdsDDD1 pdsDDD2 provides various diagonal matrices
% for each value of LSproblem.

%-----------------------------------------------------------------------
% 17 Mar 1998: First version to go with pdsco.m and lsqr.m.
% 01 Apr 1998: global pdsDDD1 pdsDDD2 now used for efficiency.
% 11 Feb 2000: Added diagonal preconditioning for LSQR, LSproblem = 1.
% 14 Dec 2000: Added diagonal preconditioning for LSQR, LSproblem = 12.
% 30 Jan 2001: Added diagonal preconditioning for LSQR, LSproblem = 21.
% 12 Feb 2001: Included in pdsco.m as private function.
%              Specialized to allow only LSproblem = 1.
%-----------------------------------------------------------------------

global pdsDDD1 pdsDDD2;

explicit  = rw(1);
LSproblem = rw(2);
LSmethod  = rw(3);
LSdamp    = rw(4);
precon    = rw(7);

if LSproblem == 1
    % The operator M is [ DA';  LSdamp*I ].
    m = nlsqr;
    n = mlsqr - m;
    if mode == 1
       if precon, x = pdsDDD2 .* x; end
       %t = feval( Aname, 2, m, n, x );   % Ask 'aprod' to form t = A'x.
       t = pdsxxxmat( Aname, 2, m, n, x );   % Ask 'aprod' to form t = A'x.
       %t = Tlsqrmat( 2, m, n, x, Aname );   % Ask 'aprod' to form t = A'x.
       y = [ (pdsDDD1.*t); (LSdamp*x) ];
    else
       t = pdsDDD1.*x(1:n);
       %y = feval( Aname, 1, m, n, t );   % Ask 'aprod' to form y = At.
       y = pdsxxxmat( Aname, 1, m, n, t );   % Ask 'aprod' to form y = A t.
       %y = Tlsqrmat( 1, m, n, t, Aname );   % Ask 'aprod' to form y = A t.
       y = y  +  LSdamp * x(n+1:mlsqr);
       if precon, y = pdsDDD2 .* y; end
    end
else
    %disp('Error in pdsxxx: Only LSproblem = 1 is allowed')
    %keyboard
    error('Error in pdsxxx: Only LSproblem = 1 is allowed')
end

% 25 Jan 2003: HKH - Avoid using global matrix - Use 1st input Aname
%-----------------------------------------------------------------------
% End private function pdsxxxlsqrmat
%-----------------------------------------------------------------------


function y = pdsxxxmat( Aname, mode, m, n, x )
%        y = pdsxxxmat( Aname, mode, m, n, x )
%    computes y = Ax (mode=1) or A'x (mode=2)
%    for some matrix A, for use with pdsco.m.

%-----------------------------------------------------------------------
% 04 Apr 1998: Default A*x and A'*y function for pdsco.m.
%              It assumes A is the global matrix pdsAAA created by pdsco.m
%              from the user's input parameter A.
% 25 Jan 2003: HKH - Avoid using global matrix - Use 1st input Aname
%-----------------------------------------------------------------------

% Original:§ function y = pdsxxxmat( mode, m, n, x )

if ischar(Aname)
   y = feval( Aname, mode, m, n, x );
else
   if mode==1,  y = Aname*x;  else  y = Aname'*x;  end
end

% Original pdsco code
% global pdsAAA;
% if mode == 1
%    y = pdsAAA*x;
% else
%    y = pdsAAA'*x;
% end

%-----------------------------------------------------------------------
% End private function pdsxxxmat
%-----------------------------------------------------------------------

function [ x, istop, itn, outfo ] = ...
   pdsTlsqr( m, n, aprodname, iw, rw, b, damp, ...
              atol, btol, conlim, itnlim, show, info )

% Call Tomlab LSQR, Tlsqr. Use special PDCO/PDSCO interface to Tlsqr
% atol algorithm implemented in Tlsqr
% No callbacks if iw has A matrix, pdsxxxlsqrmat implemented in Tlsqr
% If iw string, use pdTlsqrmat to call feval(iw, ...)

% Standard call to Tlsqr:
%
% function [ x, iStop, Iter, rNorm, xNorm, StdErr, aNorm, aCond, arNorm ] =  ...
%     Tlsqr( m, n, Aname, iw, rw, b, damp, aTol, bTol, condLim, MaxIter, ...
%            WantStdErr, nOut, D, x_0 )
%
% Special PDCO interface has the following parameters
%
% function [ x, iStop, Iter, atolold, atolnew, r3ratio ] =  ...
%     Tlsqr( m, n, Aname, iw, rw, b, damp, aTol, bTol, condLim, MaxIter, ...
%            WantStdErr, nOut, D, x_0, r3norm, atolmin, D1, D2 )
%

global pdsDDD1 pdsDDD2;

if ischar(iw)
   [ x, istop, itn, atolold, atolnew, r3ratio ] =  ...
        Tlsqr( m, n, aprodname, iw, rw, b, damp, ...
        atol, btol, conlim, itnlim, show, 0, [], [], info.r3norm, info.atolmin);
else
   if rw(7) == 1 % Precon true
      [ x, istop, itn, atolold, atolnew, r3ratio ] =  ...
           Tlsqr( m, n, iw, iw, rw, b, damp, ...
           atol, btol, conlim, itnlim, show, 0, pdsDDD2, [], ...
           info.r3norm, info.atolmin, pdsDDD1, rw(4));
   else
      [ x, istop, itn, atolold, atolnew, r3ratio ] =  ...
           Tlsqr( m, n, iw, iw, rw, b, damp, ...
           atol, btol, conlim, itnlim, show, 0, [], [], ...
           info.r3norm, info.atolmin, pdsDDD1, rw(4));
   end
end

outfo.atolold = atolold;
outfo.atolnew = atolnew;
outfo.r3ratio = r3ratio;
%-----------------------------------------------------------------------
% End private function pdTlsqr
%-----------------------------------------------------------------------
