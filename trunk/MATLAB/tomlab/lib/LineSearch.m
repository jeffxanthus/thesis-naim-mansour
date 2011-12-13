% Modified version of Fletchers line search algorithm in
% Roger Fletcher: Practical Optimization, 2nd ed., 1987.
%
% function Result = LineSearch( f, g, x, p, f_0, g_0, LineParam, alphaMax,...
%                               pType, PriLev, Prob, varargin);
%
% INPUT:
%
% f     Name of function value routine
% g     Name of gradient routine
% x     Current iterate x
% p     Search step p
% f_0   Function value at alpha=0
% g_0   Gradient at alpha=0, the directed derivative at the present point
%
% LineParam Structure. Fields used in LineParam:
%
%   LineAlg LineAlg=0 Quadratic interpolation.
%           LineAlg=1 Cubic interpolation.
%           LineAlg=2 Curvilinear Quadratic interpolation.
%           LineAlg=3 Curvilinear Cubic interpolation
%
%
%   InitStepLength Initial step length
%   sigma   rho < sigma < 1. Convergence tol. Determines accuracy
%           sigma = 0.9 inexact line search. sigma = 0.1 exact line search
%   fLowBnd    Lower bound on function at optimum
%   rho     rho value, determines rho-line, default 0.01
%   tau1    tau(1) value, how fast step grows in phase 1, default 9
%   tau2    tau(2) value, how near end point of [a,b], default 0.1
%   tau3    tau(3) value, choice in [a,b] phase 2, default 0.5
%   eps1    Minimal length for interval [a,b], normal value 1E-7
%   eps2    Minimal reduction, default 1E-12
%--------
% alphaMax  Maximal value of step length alpha
% pType     0= Normal problem. f_alpha, g_alpha defined
%           1= Nonlinear least squares. Also fields r_k = residual vector r and
%              J_k = the Jacobian matrix J are defined.  f=0.5*r'*r. g=r'*J.
%           2= Constrained nonlinear least squares.
%              Fields r_k, J_k, c_k, cJac are defined
%           3= Merit function minimization. Fields f_k, c_k, g_k, cJac defined
%           4= Penalty function minimization.
%              Defined fields f_k, c_k, g_k, cJac defined
% PriLev    If PriLev > 0 write a huge of output in LineSearch
%           If PriLev > 3 write a huge amount of output in intpol2 and intpol3
% Prob      Problem structure
%
% Extra input arguments are sent to the function and gradient routines
%
% OUTPUT:
% Result. Structure. Fields used:
%   alpha     Optimal line search step alpha.
%   f_alpha   Optimal function value at line search step alpha.
%   g_alpha   Optimal gradient value at line search step alpha.
%   alphaVec  Vector of trial step length values
%   r_k       Residual vector if NNLS problem, otherwise empty.
%   J_k       Jacobian matrix if NNLS problem, otherwise empty.
%   f_k       Function value at x+alpha*p
%   g_k       Gradient value at x+alpha*p
%   c_k       Constraint value at x+alpha*p
%   cJac      Constraint gradient value at x+alpha*p

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1994-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written April 17, 1994. Last modified Aug 13, 2009.

function Result = LineSearch( f, g, x, p, f_0, g_0, LineParam, alphaMax,...
                              pType, PriLev, Prob, varargin)

global xLast                          
% Initialization of Result structure
Result.alpha   = 0;   Result.alphaVec = 0;
Result.f_alpha = f_0; Result.g_alpha  = g_0;
Result.f_k     = [];  Result.g_k      = [];
Result.c_k     = [];  Result.cJac     = []; Result.Ax = [];

p=p(:);
x=x(:);
if nargin < 10
   PriLev=[];
   if nargin < 9
      pType=[];
      if nargin < 8
         alphaMax=[];
         if nargin < 7
            LineParam=[];
         end
      end
   end
end

if isempty(PriLev),  PriLev=0;    end
if isempty(pType),   pType=0;        end
if isempty(alphaMax),alphaMax=10; end

if isempty(LineParam)
   LineParam = LineParamDef;
end

Result.alpha = 0;

if isfield(LineParam,'r_k')
   Result.r_k     = LineParam.r_k;
else
   Result.r_k     = [];
end
if isfield(LineParam,'J_k')
   Result.J_k      = LineParam.J_k; 
else
   Result.J_k      = []; 
end

[fType gType] = CheckFunc( f, g, pType);

fLowBnd = LineParam.fLowBnd; % Lower bound on optimal function value
if isempty(fLowBnd), fLowBnd = -realmax; end
sigma  = DefPar(LineParam,'sigma',0.9); % sigma = 0.9 is inexact line search 

rho  = LineParam.rho;     % Deterines rho-line
tau1 = LineParam.tau1;    % How fast the step length in phase 1 grows
tau2 = LineParam.tau2;    % How near the end points of the interval
tau3 = LineParam.tau3;    % Choice in [a b]. (phase 2)
eps1 = LineParam.eps1;    % minimal length for the interval [a b]
eps2 = LineParam.eps2;    % minimal reduction

% Initial step
alpha_1 = LineParam.InitStepLength;
LineAlg = LineParam.LineAlg;        % Default use quadratic interpolation
MaxIter = LineParam.MaxIter;        % Maximal number of linesearch iterations

linalg=LineAlg-2*(LineAlg >=2);

%   CURV -   0 = normal line search. 1 = Curvilinear search
CURV=LineAlg >=2;    % 0 = normal line search. 1 = Curvilinear search

if CURV~=0
   if any(xLast==Inf)
      CURV=0;
   end
end

%	Initialization

k = 1;	% iteration counter
if isempty(f_0) 
   Prob.Mode = 2; % Assume also g_0 is unknown now
   [f_0, Result] = evalFunc( f, x, fType, Result, Prob, varargin{:});
   Result.f_alpha = f_0;
end

check=0;
if isempty(g_0)
   Prob.Mode = 1;
   [g_0, Result] = evalGrad( g, x, gType, Result, Prob, varargin{:});
   check=1;
   Result.g_alpha = g_0;
end

alphaVec=[];
f_alpha_km1 = f_0;
g_alpha_km1 = g_0;
fp_0 = g_0'*p;			% derivative in x in search direction,
                                % fp(x) =  g(x)' * p

if check
   if fp_0 >=0 % No descent
      Result.alpha=-100;
      Result.alphaVec=alphaVec;
      Result.f_alpha=f_0;
      Result.g_alpha=g_0;
      return; % Default values in Result are returned
   end
end
if fp_0 == 0
   my = alphaMax; % maximal value for alpha
   if PriLev > 0
      fprintf('ERROR IN LINE SEARCH! Directed derivative == 0. ');
      fprintf('alphaMax = %14.12f\n',alphaMax);
   end
else
   my = min(alphaMax,(fLowBnd - f_0) / (rho*fp_0)); % maximal value for alpha
end
if my < 0
   my=alphaMax;
   if PriLev > 0
      fprintf('ERROR IN LINE SEARCH! Upper bound < 0. Change to %15.12f\n',my);
      fprintf('Directed Derivative (should be < 0 ): %14.12f\n',fp_0);
   end
   if fp_0 > 0
      fp_0 =-fp_0;
   end
end
fp_alpha_km1 = fp_0;
alpha_k = min([1 my alpha_1]);	% alpha k   (alpha_1)

alpha_km1 = 0;			% alpha k-1 (alpha = 0)
if alpha_k < 1E-14
   if PriLev > 0
      fprintf('LineSearch: Too small step %16.12f. Set alpha==0\n',alpha_k)
   end
   Result.alpha=0;
   Result.alphaVec=alphaVec; % Either return empty or 0
   Result.f_alpha=f_0;
   Result.g_alpha=g_0;
   return  % Default values in Result are returned
end

Result1=Result;
Result=[];

%------------------------------------------------------------------------------
%
%	PHASE 1	- "bracketing phase"	(alpha  called a  in PHASE 1)
%                                            k         k
%if PriLev > 0
%	disp('alpha-max my :');
%	disp(my);
%end
while 1
   if CURV==0
      xEval=x+alpha_k*p;
   else
      xEval=x+alpha_k*p+alpha_k^2*(xLast-x+p);
   end

   Prob.Mode = 0;
   [f_alpha_k, Result] = evalFunc( f, xEval, fType, Result,Prob, varargin{:});

   while alpha_k > 1E-14 & (isinf(f_alpha_k) | isnan(f_alpha_k))
      % Try to get feasible in the number system
      alpha_k=alpha_k*0.01;
      if CURV==0
         xEval=x+alpha_k*p;
      else
         xEval=x+alpha_k*p+alpha_k^2*(xLast-x+p);
      end
      Prob.Mode = 0;
      [f_alpha_k, Result] = evalFunc( f,xEval, fType, Result,Prob,varargin{:});
      if f_alpha_k > 1E50*f_0
         % Try reducing more
         %alpha_k=alpha_k*10;
         f_alpha_k=Inf;
      end
   end
   alphaVec=[alphaVec, alpha_k];

% End of line search if f(a ) is less then fLowBnd
%                          k
   if f_alpha_k <= fLowBnd
      if PriLev > 0
         disp('--- End of line search, f less or equal than fLowBnd ---');
      end
      Prob.Mode = 1;
      [g_alpha_k, Result] = evalGrad(g, xEval, gType, Result,Prob, varargin{:});

      Result.alpha=alpha_k;
      Result.alphaVec=alphaVec;
      Result.f_alpha=f_alpha_k;
      Result.g_alpha=g_alpha_k;
      return;			% End of line search
   end
%
% End of PHASE 1 if f(a ) > f(a   )  or f(a ) goes over the rho-line
%                      k       k-1         k
   if f_alpha_k >=f_alpha_km1 | ...
      f_alpha_k > f_0+alpha_k*rho*fp_0 % Test 2.5.1 (rho-line)

      if f_alpha_k == f_alpha_km1 & norm(p) <= 1E-8 % Extremely small step
         if PriLev > 0
           disp('! End of line search, f(alpha)==f(alpha_km1). |p| < 1E-8');
         end
         Prob.Mode = 1;
         [g_alpha_k, Result] = evalGrad(g,xEval,gType,Result,Prob,varargin{:});

         Result.alpha=alpha_k;
         Result.alphaVec=alphaVec;
         Result.f_alpha=f_alpha_k;
         Result.g_alpha=g_alpha_k;
         return;
      end
      a = alpha_km1;		% Interval [a b] for PHASE 2
      b = alpha_k;
      f_a = f_alpha_km1;
      g_a = g_alpha_km1;
      fp_a = fp_alpha_km1;
      f_b = f_alpha_k;
      R=Result1;
      Result1=Result;
      Result=R;
      if linalg > 0
% Evaluate f'(b = new worse point). Only needed because of cubic interpolation.
         Prob.Mode = 1;
         [g_b, Result1] = evalGrad(g, xEval, gType, Result1,Prob, varargin{:}); 
         fp_b = g_b'*p;
      end
      if PriLev > 0
disp('--- End PHASE 1, f(alpha) > rho-line or f(alpha) >= prev.f(alpha). ---');
      end
      break;
   end
% Evaluate f'(alpha_k)
   Prob.Mode = 1;
   [g_alpha_k, Result] = evalGrad(g, xEval, gType, Result,Prob, varargin{:});

   fp_alpha_k = g_alpha_k'*p;
%
% End line search if Fletchers 2nd condition is fulfilled
%
   if abs(fp_alpha_k) <= -sigma*fp_0 | ...
      (fp_alpha_k < 0 & alpha_k==my)
      if PriLev > 0
         fprintf('f_alpha_k = %16.12f alpha=%16.12f\n',f_alpha_k,alpha_k);
         fprintf('fp_alpha_k = %16.12f fp_0=%16.12f\n',fp_alpha_k,fp_0);
         disp('--- End of line search, F-2 in Phase I  ---');
      end
      Result.alpha=alpha_k;
      Result.alphaVec=alphaVec;
      Result.f_alpha=f_alpha_k;
      Result.g_alpha=g_alpha_k;
      return;
   end
%
% End of PHASE 1 if the function grows in the search direction
%
   if fp_alpha_k >= 0
      a = alpha_k;		% Turn the interval [a b]
      b = alpha_km1;
      f_a = f_alpha_k;
      f_b = f_alpha_km1;
      g_a = g_alpha_k;
      g_b = g_alpha_km1;
      fp_a = fp_alpha_k;
      fp_b = fp_alpha_km1;
      if PriLev > 0
         disp('--- End of phase 1, f grows in the search direction ---');
      end
      break;
   end
%
% Update of a    and a
%            k-1      k
   lower = 2*alpha_k - alpha_km1;
   if my <= lower;			% Determine next alpha
      alpha_kp1 = my;			% Maximal alpha
   else
      upper = min(my, alpha_k+tau1*(alpha_k-alpha_km1));
					% Choose alpha(k+1) in [lower, upper]
      if upper < 1E-14
         if PriLev > 0
            fprintf('LineSearch: upper limit too low in PHASE 1, %d its.',k)
            fprintf('upper %16.8e. Set alpha==0\n',upper)
         end
         Result.alpha=alpha_k;
         Result.alphaVec=alphaVec;
         Result.f_alpha=f_alpha_k;
         Result.g_alpha=g_alpha_k;
         return % Default values in Result are OK
      end
      if linalg == 0
         alpha_kp1 = intpol2(alpha_k, f_alpha_k, fp_alpha_k,   ...
                          alpha_km1, f_alpha_km1, lower, upper, PriLev);
      elseif linalg==1
         alpha_kp1 = intpol3(alpha_k, f_alpha_k, fp_alpha_k,  alpha_km1, ...
                           f_alpha_km1, fp_alpha_km1, lower, upper, PriLev);
      else
         alpha_kp1 = intpol2(alpha_k, f_alpha_k, fp_alpha_k,   ...
                          alpha_km1, f_alpha_km1, lower, upper, PriLev);
      end
   end	% (if)
   alpha_km1 = alpha_k;
   f_alpha_km1 = f_alpha_k;
   g_alpha_km1 = g_alpha_k;
   fp_alpha_km1 = fp_alpha_k;
   alpha_k = alpha_kp1;
   k = k + 1;

   if k > MaxIter 
      if PriLev > 0
         fprintf('alpha_k: %f, f_alpha_k %f\n',alpha_k,f_alpha_k)
         fprintf('LineSearch: Too many PHASE 1 iterations, %d its.\n',k)
      end
      if isempty(Result), Result=Result1; end
      Result.alpha=alpha_km1;    % Not allowed to compute f for new alpha_k
      Result.alphaVec=alphaVec;
      Result.f_alpha=f_alpha_km1;
      Result.g_alpha=g_alpha_km1;
      return % Default values in Result are OK
   end
   Result1=Result;
   Result=[];

   if PriLev > 0
      if k == 1
         fprintf('k = %4.0f  alpha = %14.12f \n',k,alpha_k);
      end
   end
end


%------------------------------------------------------------------------------
%
%	PHASE 2	- "sectioning phase"
%
if PriLev > 0
   fprintf('PHASE 2 === Interval [a, b] =%16.12f%16.12f\n',a,b);
end

Result2=Result1;
Result1=Result;
Result=[];

while 1
%
% Extra test 1, end if the interval is too small
%
   if abs(a-b) <= eps1
      if f_a < f_b | (f_a==f_b & a > b)
         alpha_k = a;
         f_alpha_k = f_a;
         g_alpha_k = g_a;
         Result=Result1;
      else
         Result=Result2;
         alpha_k = b;
         f_alpha_k = f_b;
         if exist('g_b','var')
            g_alpha_k = g_b;
            if isfield(Result2,'J_k')
               Result.J_k = Result2.J_k; % HKH 001110 This line seems needed
            end
         else                
            % Compute g_b if it hasn't been computed before
            Prob.Mode = 1;
            [g_b, Result] = evalGrad(g, xEval, gType, Result,Prob,varargin{:});
            g_alpha_k = g_b;
         end
      end
      if PriLev > 0
         fprintf('=== End line search, [a b] < eps1 = %16.12f.',eps1);
         fprintf(' alpha_k = %16.12f\n',alpha_k);
      end
      Result.alpha=alpha_k;
      Result.alphaVec=alphaVec;
      Result.f_alpha=f_alpha_k;
      Result.g_alpha=g_alpha_k;
      return;
   end
%
% Compute new alpha using interpolation
%
   lower = a + tau2*(b-a);
   upper = b - tau3*(b-a);		% Choose alpha in [lower upper]
   if linalg == 0
      alpha_k = intpol2(a, f_a, fp_a, b, f_b, lower, upper, PriLev);
   elseif linalg==1
      alpha_k = intpol3(a, f_a, fp_a, b, f_b, fp_b, lower, upper, PriLev);
   else
      alpha_k = intpol2(a, f_a, fp_a, b, f_b, lower, upper, PriLev);
   end	% (if)
   if PriLev > 0
      fprintf('Interpolated alpha, alpha_k = %16.12f\n',alpha_k);
   end
%
% Extra test 2. Check if reduction possible.
%
   if (a-alpha_k)*fp_a <= eps2
      if PriLev > 0
         fprintf('alpha_k=%18.12f; a=%18.12f',alpha_k,a);
         fprintf(' END: Too small descent\n');
      end
      Result=Result1;
      alpha_k = a;
      f_alpha_k = f_a;
      g_alpha_k = g_a;
      Result.alpha=alpha_k;
      Result.alphaVec=alphaVec;
      Result.f_alpha=f_alpha_k;
      Result.g_alpha=g_alpha_k;
      return;
   end
%
% Compute new f(alpha)
%
   if CURV==0
      xEval=x+alpha_k*p;
   else
      xEval=x+alpha_k*p+alpha_k^2*(xLast-x+p);
   end
   Prob.Mode = 2*(linalg > 0); % signal if we know of derivative call
   [f_alpha_k, Result] = evalFunc( f, xEval, fType, Result,Prob, varargin{:});

   alphaVec=[alphaVec, alpha_k];
   if PriLev > 0
      fprintf('OLD BEST f(a) = %16.12f    a=%16.12f \n',f_a,a);
      fprintf('New  f_alpha_k = %16.12f alpha=%16.12f\n',f_alpha_k,alpha_k);
   end
%
% New alpha better?
%
   if f_alpha_k >= f_a | f_alpha_k > f_0 + rho*alpha_k*fp_0
      b = alpha_k;		% New alpha_k worse, move b
      f_b = f_alpha_k;		% 
      Result2=Result;

% Evaluate f'(b) only because of cubic interpolation!!
      if linalg > 0
         Prob.Mode = 1;
         [g_b, Result2] = evalGrad(g, xEval, gType, Result2,Prob, varargin{:});

         fp_b = g_b'*p;
      end
   else				% New alpha_k better than a
      Prob.Mode = 1;
      [g_alpha_k, Result] = evalGrad(g, xEval, gType, Result,Prob, varargin{:});

      fp_alpha_k = g_alpha_k'*p;
      if abs(fp_alpha_k) <= -sigma*fp_0
         if PriLev > 0
            fprintf('f_alpha_k = %16.12f alpha=%16.12f\n',f_alpha_k,alpha_k);
            fprintf('fp_alpha_k = %16.12f fp_0=%16.12f\n',fp_alpha_k,fp_0);
            disp('--- End of line search, F-2 ---');
         end
         Result.alpha=alpha_k;
         Result.alphaVec=alphaVec;
         Result.f_alpha=f_alpha_k;
         Result.g_alpha=g_alpha_k;
         return;
      end	% (if)
%
% Determine end point that shall be deleted
%
      if (b-a)*fp_alpha_k >= 0	% a is new b
         b = a;
         f_b = f_a;
         g_b = g_a;
         fp_b = fp_a;
         Result2=Result1;
      end	% (if)
      a = alpha_k;		% New alpha the best uptil now
      f_a = f_alpha_k;
      g_a = g_alpha_k;
      fp_a = fp_alpha_k;
      Result1=Result;
      Result=[];
   end	% (else)
end	% (while)

% =====================================================================
% function  alfa = intpol2(x0, f0, g0, x1, f1, a, b, PriLev)
% =====================================================================
%
% Quadratic interpolation of the function f using functions values f0 & f1 and
% the derivative g0 in the points x0. Returns the minimum of the
% interpolated second degree polynomial p in the interval [a,b] (alt. [b,a]).
%
% If PriLev > 3 write a lot of output
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Written May 22, 1995.  Last modified June 21, 1999.
%

% q(z) = f0 +  g0*z + (f1 - f0 - g0)z^2;  z in [0,1]
%
% Use alfa = a + z*(b-a) and df/dz = df/d_alfa * d_alfa/dz = df/d_alfa * (b-a)
% Then g0 = g0 * (b-a)
%
% dq/dz = g0 + 2 ( f1 - f0 - g0) z, gives dq/dz = 0 when 
%     z = -g0 / (2(f1 - f0 - g0))
% Now [x0, x1] are a bit bigger than the given [a,b], due to safe guarding.

function  alfa = intpol2(x0, f0, g0, x1, f1, a, b, PriLev)

if  nargin < 8, PriLev=0;end

if PriLev > 3
   fprintf('Quad input x0=%15.12f, a=%14.12f, b=%14.12f x1=%14.12f\n',x0,a,b,x1);
   fprintf('Quad func  f0 %15.12f,f1=%14.12f,g0=%14.12f\n',f0,f1,g0);
end
if isinf(f0)
   error('intpol2: Gradient f0 has infinite value')
end
if isinf(g0)
   error('intpol2: Gradient g0 has infinite value')
end
% Transform derivative to [0,1]
g0 = g0 * (x1-x0);
% Let c = (f1 - f0 - g0)
c = f1 - f0 - g0;
% Compute minimum
if c==0
   if PriLev > 3
      fprintf('intpol2: ERROR - computing quadratic interpolation\n')
      fprintf('g0 = %15.12f, c  = %15.12f ',g0,c)
      fprintf('f1 = %15.12f, f0 = %15.12f \n',f1,f0)   
   end
   alfa=min(a,b);
   return
else
   z = (-g0/c)/2;
end
% Function value at minmium
q_z = f0 + g0*z + c*z^2;
% Compute function value at endpoints in interval [a,b]
% Must transform [a,b] to [0,1] ( = [x0,x1] )
% z_new = (alfa-a)/(b-a) = (?-x0)/x1-x0)
z_a = (a-x0)/(x1-x0);
q_a = f0 + g0*z_a + c*z_a^2;
z_b = (b-x0)/(x1-x0);
q_b = f0 + g0*z_b + c*z_b^2;
if PriLev > 3
   fprintf('Quad minimum  %14.12f, a=%15.12f, b=%15.12f\n',z,z_a,z_b);
   fprintf('Quad func.val %14.12f,fa=%15.12f,fb=%15.12f\n',q_z,q_a,q_b);
end
% Check if z in [a,b]
inside=1;
if z < min(z_a,z_b) | z > max(z_a,z_b)
   inside=0;
   if PriLev > 3
      disp('Minimum is outside range [a,b]')
   end
end
if q_z <= q_a & inside
   if q_z <= q_b
      alfa=x0+z*(x1-x0);
   else
      alfa=x0+z_b*(x1-x0);
   end
elseif q_a < q_b
   alfa=x0+z_a*(x1-x0);
else
   alfa=x0+z_b*(x1-x0);
end
if PriLev > 3
   fprintf('Best alfa %14.12f\n',alfa);
end

% MODIFICATION LOG intpol2:
%
% 980717   mbk      'alfa=min(a,b)' and 'return' on line ~= 41 moved outside
%                   'PriLev' condition. ;
% 990615   hkh      Change exit to return in case of inf values
%

% =====================================================================
% function  alfa = intpol3(x0, f0, g0, x1, f1, g1 , a, b, PriLev)
% =====================================================================
%
% Cubic interpolation of the function f using functions values f0 & f1 and
% the derivatives g0 & g1 in the points x0 and x1. Returns the minimum of the
% interpolated third degree polynomial p in the interval [a,b] (alt. [b,a]).
%
% If PriLev > 3 write a lot of output
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Written May 22, 1995.  Last modified June 21, 1999.
%

% c(z) = f0 +  g0*z + r * z^2 + s * z^3 ;  z in [0,1]
% where r = 3*(f1-f0) - 2*g0 - g1 and s = g0 + g1 - 2*(f1 - f0)
%
% Use alfa = a + z*(b-a) and df/dz = df/d_alfa * d_alfa/dz = df/d_alfa * (b-a)
% Then g0 = g0 * (b-a) and g1 = g1 * (b-1)
%
% dq/dz = g0 + 2 r z + 3 s z^2 , gives dq/dz = 0 when 
%     z = (-r +- sqrt(r^2 -  3 s go)) / 3 s
% Now [x0, x1] are a bit bigger than the given [a,b], due to safe guarding.

function  alfa = intpol3(x0, f0, g0, x1, f1, g1 , a, b, PriLev)

if  nargin < 9, PriLev=0;end

if PriLev > 3
   fprintf('Cubic input x0=%16.12f, a=%14.12f, b=%14.12f x1=%14.12f\n',...
            x0,a,b,x1);
   fprintf('Cubic func  f0 %16.12f,f1=%14.12f,g0=%14.12f, g1=%14.12f\n',...
            f0,f1,g0,g1);
end
% Transform derivative to [0,1]
g0 = g0 * (x1-x0);
g1 = g1 * (x1-x0);
% Compute r and s
r = 3*(f1-f0) - 2*g0 - g1;
s = g0 + g1 - 2*(f1 - f0);

if abs(s) < 1E-12 * abs(r);
   % s is too small. Assume this is a quadratic instead
   if PriLev > 3
      fprintf('s, z^3-term is too small; s %14.12f r %14.12f \n',s,r);
   end
   g0 = g0/(x1-x0); % Transform derivative back to [x0,x1].
   alfa = intpol2(x0, f0, g0, x1, f1, a, b);
else
   % OK cubic term
   % Compute the two minimum points z_1 and z_2
   root = sqrt(r*r - 3 * s * g0);
   z_1 = (-r + root)/3*s;
   z_2 = (-r - root)/3*s;
   % Function value at minium
   c_1 = f0 + g0*z_1 + r*z_1^2 + s*z_1^3;
   c_2 = f0 + g0*z_2 + r*z_2^2 + s*z_2^3;
   if PriLev > 3
      fprintf('Cubic minimum 1 %14.12f, z_1= %14.12f\n',c_1,z_1);
      fprintf('Cubic minimum 2 %14.12f, z_2= %14.12f\n',c_2,z_2);
      fprintf('Min in oldCoord %14.12f,    = %14.12f\n',...
               x0+z_1*(x1-x0), x0+z_2*(x1-x0));
   end

   % Compute function value at endpoints in interval [a,b]
   % Must transform [a,b] to [0,1] ( = [x0,x1] )
   % z_new = (alfa-a)/(b-a) = (?-x0)/x1-x0)
   z_a = (a-x0)/(x1-x0);
   c_a = f0 + g0*z_a + r*z_a^2 + s*z_a^3;
   z_b = (b-x0)/(x1-x0);
   c_b = f0 + g0*z_b + r*z_b^2 + s*z_b^3;
   if PriLev > 3
      fprintf('Interval in z:  a= %14.12f, b= %14.12f\n',z_a,z_b);
      fprintf('Cubic func.val fa= %14.12f,fb= %14.12f\n',c_a,c_b);
   end
   % Check if z_1 and z_2 in [a,b] and which is best
   c_z = 1e10;
   if z_1 < min(z_a,z_b) | z_1 > max(z_a,z_b)
      inside=0;
      if PriLev > 3
         disp('Minimum z_1 is outside range [a,b]')
      end
   else
      inside=1;
      z = z_1;
      c_z = c_1;
   end
   if z_2 < min(z_a,z_b) | z_2 > max(z_a,z_b)
      if PriLev > 3
         disp('Minimum z_2 is outside range [a,b]')
      end
   else
      if inside % Check which function value is best
         if c_2 < c_z
            z = z_2;
            c_z = c_2;
            if PriLev > 3
               disp('Minimum z_2 is BEST')
            end
         else
            if PriLev > 3
               disp('Minimum z_1 is BEST')
            end
         end
      else
         inside=1;
         z = z_2;
         c_z = z_2;
      end
   end
   if c_z <= c_a & inside
      if c_z <= c_b
         alfa=x0+z*(x1-x0);
      else
         alfa=x0+z_b*(x1-x0);
      end
   elseif c_a < c_b
      alfa=x0+z_a*(x1-x0);
   else
      alfa=x0+z_b*(x1-x0);
   end
   if PriLev > 3
      fprintf('Best alfa %14.12f\n',alfa);
   end
end

% MODIFICATION LOG intpol3:
%
% 980717   mbk      Transform derivative g0 back to [x0,x1] before calling
%                   inpol2.
%                   'r' changed to '-r' when computing 'z_1' and 'z_2'.
%
% ============================================================
function [fType,gType] = CheckFunc( f, g, pType)
% ============================================================
if strcmp(f,'nlp_f')
   fType=1;
elseif strcmp(f,'nlp_r')
   fType=2;
elseif strcmp(f,'con_fm')
   fType=3;
else
   fType=4+pType;
end

if strcmp(g,'nlp_g')
   gType=1;
elseif strcmp(g,'nlp_J')
   gType=2;
elseif strcmp(g,'con_gm')
   gType=3;
else
   gType=4+pType;
end


% ============================================================
function [f_k, Result] = evalFunc( f, x, fType, Result, Prob, varargin)
% ============================================================

nargin;
switch fType
case 1
   f_k = nlp_f(x,Prob, varargin{:}); 
case 2
   r_k = nlp_r(x , Prob, varargin{:} );     % Function value
   f_k = 0.5 * (r_k' * r_k);
   Result.r_k = r_k;
case 3
   [f_k Result] = con_fm(x, Result, Prob, varargin{:}); 
   % Set in con_fm: Result.f_k Result.c_k Result.Ax
case 4
   f_k = feval(f,x,Prob, varargin{:});
case {5,6}
   r_k = feval(f,x,Prob, varargin{:});
   f_k = 0.5 * (r_k' * r_k);
   Result.r_k = r_k;
case {7,8}
   [f_k Result] = feval(f, x, Result, Prob, varargin{:}); 
   % Set in con_fm: Result.f_k Result.c_k Result.Ax
end
   
% ============================================================
function [g_k, Result] = evalGrad( g, x, gType, Result, Prob, varargin)
% ============================================================

nargin;
switch gType
case 1
   g_k = nlp_g(x, Prob, varargin{:}); 
case 2
   J_k = nlp_J(x , Prob, varargin{:} );
   g_k = J_k' * Result.r_k;
   Result.J_k = J_k;
case 3
   [g_k Result] = con_gm(x, Result, Prob, varargin{:}); 
   % Set in con_gm: Result.g_k Result.cJac
case 4
   g_k = feval(g, x, Prob, varargin{:}); 
case {5,6}
   J_k = feval(g, x, Prob, varargin{:}); 
   g_k = J_k' * Result.r_k;
   Result.J_k = J_k;
case {7,8}
   [g_k Result] = feval(g, x, Result, Prob, varargin{:}); 
   % Set in con_gm: Result.g_k Result.cJac
end

% MODIFICATION LOG:
%
% 980918  hkh  Changed f_min to fLowBnd
% 980919  hkh  Changed to use input structure optParam instead of vector param
% 981005  hkh  Change to more general structure Result as output, including
%              line search on merit functions and penalty functions.
%              Changing function name from linesrch.m to LineSearch.m.
%              Changing variable names.
% 981017  hkh  Errors in sending the correct c_k to g (pType==3,4).
%              Changing the logic for the merit function computations
% 981019  hkh  Not setting Result.f_k to f_0 at start, error if alpha=0 and
%              merit function.
% 981027  hkh  Move fLowBnd to optParam.fLowBnd
% 981110  hkh  Change dc_k to cJac in Result and OUT structure
%              Change to alpha in comments and printings
% 000212  hkh  Correct handling of Result structure. Wrong r_k sometimes.
% 000923  hkh  Change to use all parameters from structure LineParam
% 001206  hkh  Init with Result.alpha = 0.
% 010223  hkh  Safeguarding if max number of iterations in phase 1
% 011204  hkh  Some information not returned if maximal number of iterations
% 020409  hkh  Set Prob.Mode, send Prob structure explicit
% 041018  hkh  Allow LineParam.sigma [], set to 0.9 as default
% 090813  med  mlint check