% eloLine
%
% Line search algorithm for efficient local optimization, i.e.
% when function values are costly and gradients must be estimated
% numerically
%
% function Result = eloLine( f, g, x, p, f_0, g_0, r_0, J_0, LineParam, ...
%                            alphaMax, alpha_1, pType, PriLev, varargin);
%
% INPUT:
%
% f     Name of function value routine
% g     Name of gradient routine
% x     Current iterate x
% p     Search step p
% f_0   Function value at alpha=0
% g_0   Gradient at alpha=0, the directed derivative at the present point
% r_0   Residual vector at alpha=0 if NNLS problem, otherwise empty.
% J_0   Jacobian matrix at alpha=0 if NNLS problem, otherwise empty.
%
% LineParam Fields used:
%   LineAlg LineAlg=0 Quadratic interpolation.
%           LineAlg=1 Cubic interpolation.
%           LineAlg=2 Curvilinear Quadratic interpolation.
%           LineAlg=3 Curvilinear Cubic interpolation
%
%   sigma   rho < sigma < 1. Convergence tol. Determines accuracy
%           sigma = 0.9 inexact line search. sigma = 0.1 exact line search
%   fLow    Lower bound on function at optimum
%   rho     rho value, determines rho-line, default 0.01
%   tau1    tau(1) value, how fast step grows in phase 1, default 9
%   tau2    tau(2) value, how near end point of [a,b], default 0.1
%   tau3    tau(3) value, choice in [a,b] phase 2, default 0.5
%   eps1    Minimal length for interval [a,b], normal value 1E-7
%   eps2    Minimal reduction, default 1E-12
%--------
% alphaMax  Maximal value of step length alpha
% alpha_1   1st step in alpha.
% pType     0= Normal problem. f_alpha, g_alpha defined
%           1= Nonlinear least squares. Also fields r_k = residual vector r and
%              J_k = the Jacobian matrix J are defined.  f=0.5*r'*r. g=r'*J.
%           2= Constrained nonlinear least squares.
%              Fields r_k, J_k, c_k, cJac are defined
%           3= Merit function minimization. Fields f_k, c_k, g_k, cJac defined
%           4= Penalty function minimization.
%              Defined fields f_k, c_k, g_k, cJac defined
% PriLev    If PriLev > 0 write a huge of output in eloLine
%           If PriLev > 3 write a huge amount of output in intpol2 and intpol3
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

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., Sweden. $Release: 7.2.0$
% Written Feb 5, 2000. Last modified Jul 17, 2009.

function Result = eloLine( f, g, x, p, f_0, g_0, r_0, J_0, LineParam, ...
                           alphaMax, alpha_1, pType, PriLev, varargin)

% Initialization of Result structure
Result.alpha   = 0;   Result.alphaVec = 0;
Result.f_alpha = f_0; Result.g_alpha  = g_0;
Result.r_k     = [];  Result.J_k      = [];
Result.f_k     = [];  Result.g_k      = [];
Result.c_k     = [];  Result.cJac     = []; Result.Ax = [];

p=p(:);
x=x(:);
if nargin < 11
   PriLev=[];
   if nargin < 10
      pType=[];
      if nargin < 9
         alpha_1=[];
         if nargin < 8
            alphaMax=[];
            if nargin < 7
               LineParam=[];
            end
         end
      end
   end
end

if isempty(PriLev),  PriLev=0;    end
if isempty(pType),   pType=0;        end
if isempty(alpha_1), alpha_1=1;   end
if isempty(alphaMax),alphaMax=10; end

if isempty(LineParam), LineParam=LineParamDef; end

disp('ELOLINE')
PriLev=1;

[fType gType] = CheckFunc( f, g, pType);

if ~isempty(r_0)
   NLLS=1;
   v1 = J_0*p;
else
   NLLS=0;
end

fLow = LineParam.fLowBnd; % Lower bound on optimal func value
sigma  = DefPar(LineParam,'sigma',0.9); % sigma = 0.9 is inexact line search
rho  = LineParam.rho;     % Determines rho-line
tau1 = LineParam.tau1;    % How fast the step length in phase 1 grows
tau2 = LineParam.tau2;    % How near the end points of the interval
tau3 = LineParam.tau3;    % Choice in [a b]. (phase 2)
eps1 = LineParam.eps1;    % minimal length for the interval [a b]
eps2 = LineParam.eps2;    % minimal reduction

MaxIter = LineParam.MaxIter;  % Max number of linesearch iterations

LineAlg=LineParam.LineAlg;    % Default use cubic interpolation

%   eps2,eps1,tau3,tau2,tau1,rho
%   alpha_1,fLow,sigma,alphaMax

%
%	Initialization
%
alpha_1 = abs(alpha_1);

k = 1;				% iteration counter
alpha_k=0;
if isempty(f_0)
   [f_0, Result] = evalFunc( f, x, fType, Result, varargin{:});
   Result.f_alpha = f_0;
end

check=0;
if isempty(g_0)
   [g_0, Result] = evalGrad( g, x, gType, Result, varargin{:});
   check=1;
   Result.g_alpha = g_0;
end

alphaVec=[];
f_alpha_km1 = f_0;
g_alpha_km1 = g_0;
fff=0.5*(r_0'*r_0);
if fff~=f_0
   disp('ERROR!!!!')
   f_0
   fff
   pause
end
fp_0 = g_0'*p;			% derivative in x in search direction,
                                % fp(x) =  g(x)' * p
fprintf('fp_0 %40.20f\n',fp_0);
pause
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
   end
   fprintf('alphaMax = %16.8e\n',alphaMax);
else
   my = min(alphaMax,(fLow - f_0) / (rho*fp_0)); % maximal value for alpha
end
if my < 0
   my=alphaMax;
   if PriLev > 0
      fprintf('ERROR IN LINE SEARCH! Upper bound < 0. Change to %15.3f\n',my);
      fprintf('Directed Derivative (should be < 0 ): %16.8e \n',fp_0);
   end
   if fp_0 > 0
      fp_0 =-fp_0;
   end
end
fp_alpha_km1 = fp_0;
alpha_k = min([1 my alpha_1]);	% alpha k   (alpha_1)
fp_alpha_k=Inf;

alpha_km1 = 0;			% alpha k-1 (alpha = 0)
if alpha_k < 1E-14
   if PriLev > 0
      fprintf('eloLine: Too small step %16.8e. Set alpha==0\n',alpha_k)
   end
   return  % Default values in Result are returned
end

%------------------------------------------------------------------------------
%
%	PHASE 1	- "Aggresive bracketing phase"	(alpha  called a  in PHASE 1)
%                                                     k         k

Result1=Result;
Result=[];

done = 1;		% IN PHASE 1 when done == 1;
while 1
   xEval=x+alpha_k*p;

   [f_alpha_k, Result] = evalFunc( f, xEval, fType, Result,varargin{:});

   while alpha_k > 1E-14 & (isinf(f_alpha_k) | isnan(f_alpha_k))
      % Try to get feasible in the number system
      alpha_k=alpha_k*0.01;
      xEval=x+alpha_k*p;
      [f_alpha_k, Result] = evalFunc( f, xEval, fType, Result,varargin{:});
      if f_alpha_k > 1E50*f_0
         % Try reducing more
         %alpha_k=alpha_k*10;
         f_alpha_k=Inf;
      end
   end


   alphaVec=[alphaVec, alpha_k];

%
% End of line search if f(a ) is less then fLow
%                          k
   if f_alpha_k <= fLow
      if PriLev > -100
         disp('--- End of line search, f less or equal than fLow ---');
      end
      [g_alpha_k, Result] = evalGrad(g, xEval, gType, Result,varargin{:});

      Result.alpha=alpha_k;
      Result.alphaVec=alphaVec;
      Result.f_alpha=f_alpha_k;
      Result.g_alpha=g_alpha_k;
      return;			% End of line search
   end
%
% End of PHASE 1 if f(a ) > f(a   )  or f(a ) goes over the rho-line
%                      k       k-1         k
   if f_alpha_k <= f_0+alpha_k*fp_0
      % Concave. Try aggressive strategy
      fprintf('Concave: %14.12f <= %14.12f',f_alpha_k,f_0+alpha_k*fp_0);
      alphaNew=min(alpha_k*2,alphaMax);
alphaNew
      fNew=-Inf;
      ResultNew=[];
      while alphaNew <= alphaMax & alphaNew > alpha_k & fNew < f_alpha_k
         [fNew, ResultNew] = evalFunc( f, xEval, fType, ResultNew,varargin{:});
         fprintf('Concave: New %14.12f <= %14.12f',alphaNew,fNew);
         alphaVec=[alphaVec, alphaNew];
         if fNew < f_alpha_k
            alpha=alphaNew;
            f_alpha_k=fNew;
            Result=ResultNew;
            alphaNew=min(alpha*2,alphaMax);
            ResultNew=[];
         end
      end
      [g_alpha_k, Result] = evalGrad(g, xEval, gType, Result,varargin{:});

      Result.alpha=alpha_k;
      Result.alphaVec=alphaVec;
      Result.f_alpha=f_alpha_k;
      Result.g_alpha=g_alpha_k;
      return;
   elseif f_alpha_k >=f_alpha_km1 | ...
      f_alpha_k > f_0+alpha_k*rho*fp_0 % Test 2.5.1 (rho-line)

      if f_alpha_k == f_alpha_km1 & norm(p) <= 1E-8 % Extremely small step
         if PriLev > 0
           disp('! End of line search, f(alpha)==f(alpha_km1). |p| < 1E-8');
         end

         [g_alpha_k, Result] = evalGrad(g, xEval, gType, Result,varargin{:});

         Result.alpha=alpha_k;
         Result.alphaVec=alphaVec;
         Result.f_alpha=f_alpha_k;
         Result.g_alpha=g_alpha_k;
         return;
      end
      done = 0;		% Go to PHASE 2
      if alpha_km1==0
         % Put alpha_k as a, even if worse value than 0
         b = 0;		% Interval [a b] for PHASE 2
         a = alpha_k;
         f_a = f_alpha_k;
         f_b = f_alpha_km1;
         fp_a = fp_alpha_k;
         % Result stores alpha_k, Result1 stores b=0
      else
         a = alpha_km1;		% Interval [a b] for PHASE 2
         b = alpha_k;
         f_a = f_alpha_km1;
         f_b = f_alpha_k;
         %g_a = g_alpha_km1;
         fp_a = fp_alpha_km1;
         % Numerical estimate of fp_b
         fp_b=(f_b-f_a)/(b-a);
         R=Result1;
         Result1=Result;
         Result=R;
      end

      if PriLev > 0
disp('--- End PHASE 1, f(alpha) > rho-line or f(alpha) >= prev.f(alpha). ---');
      end
      break;
   end
   if done==0
      disp('ERROR in LINESEARCH. Should be in PHASE 2');
      pause
   end
% Evaluate f'(alpha_k)
   %[g_alpha_k, Result] = evalGrad(g, xEval, gType, Result,varargin{:});
   %fp_alpha_k = g_alpha_k'*p;
   if alpha_k~=alpha_km1 & alpha_km1~=0
      fp_alpha_k=(f_alpha_k-f_alpha_km1)/(alpha_k-alpha_km1);
   disp('compute fp_alpha in Phase I')
fp_alpha_k
pause
   else
      fp_alpha_k=Inf;
   end
%
% End line search if Fletchers 2nd condition is fulfilled
%
   if abs(fp_alpha_k) <= -sigma*fp_0 | ...
      (fp_alpha_k > 0 & alpha_km1==my)
      done = 1;
      if PriLev > 0
         fprintf('f_alpha_k = %16.8e alpha=%16.8e\n',f_alpha_k,alpha_k);
         fprintf('fp_alpha_k = %16.8e fp_0=%16.8e\n',fp_alpha_k,fp_0);
         disp('--- End of line search, F-2 in Phase I  ---');
      end
      [g_alpha_k, Result] = evalGrad(g, xEval, gType, Result,varargin{:});
      Result.alpha=alpha_k;
      Result.alphaVec=alphaVec;
      Result.f_alpha=f_alpha_k;
      Result.g_alpha=g_alpha_k;
      return;
   end
%
% End of PHASE 1 if the function grows in the search direction
%
   if ~isinf(fp_alpha_k) & ( fp_alpha_k >= 0 ...
      | (alpha_km1~=0 & (f_alpha_km1-f_alpha)*(alpha_km1-alpha_k) > 0))

      done = 0;		% Go to PHASE 2
      a = alpha_k;		% Turn the interval [a b]
      b = alpha_km1;
      f_a = f_alpha_k;
      f_b = f_alpha_km1;
      fp_a = fp_alpha_k;
      fp_b = fp_alpha_km1;
      if PriLev > 0
         disp('--- End of phase 1, f grows in the search direction ---');
      end
      break;
   end
   if done==0
      disp('ERROR in LINESEARCH. Should be in PHASE 2');
      pause
   end
%
% Update of a    and a
%            k-1      k
   if my == alpha_k
      a = alpha_k;
      b = alpha_km1;
      f_a = f_alpha_k;
      f_b = f_alpha_km1;
      fp_a = fp_alpha_k;
      fp_b = fp_alpha_km1;
      break;
   end
   lower = 2*alpha_k - alpha_km1;
   if my <= lower;			% Determine next alpha
      alpha_kp1 = my;			% Maximal alpha
   else
      upper = min(my, alpha_k+tau1*(alpha_k-alpha_km1));
					% Choose alpha(k+1) in [lower, upper]
      if max(lower,upper) < 1E-11
         if PriLev > 0
            fprintf('eloLine: upper limit too low in PHASE 1, %d its.',k)
            fprintf('upper %16.8e. Set alpha==0\n',upper)
         end
         return % Default values in Result are OK
      end
      if NLLS
         % v1 = J_0*p; is computed
         %FDIFF = Result.r_k-r_0;
         %v2 = FDIFF-v1;
         FDIFF = (Result.r_k-r_0)/alpha_k;
         %xprinte(FDIFF,'FDIFF:')
         v2 = (FDIFF-v1)/alpha_k;
         [alpha_kp1] = parrm(r_0, v1, v2,alpha_k);
alpha_k
alpha_kp1
my
         fprintf('parrm-1: %16.14f in%17.14f%17.14f\n',alpha_kp1,lower,upper);
         alpha_kp1=min(max(min(lower,upper),alpha_kp1),max(upper,lower));
pause
      elseif LineAlg==0
         % Always use f'(0), and two function values, to avoid computing f'(x)
         alpha_kp1 = intpol2(0, f_0, fp_0,   ...
                          alpha_k, f_alpha_k, lower, upper, PriLev);
         fprintf('intpol2: %16.14f in%17.14f%17.14f\n',alpha,lower,upper);
      elseif LineAlg==1
         % Use f'(0) and three function values f'(0),f(a),f(b)
         % MUST REWRITE INTPOL3, fp_b to a,f_a
         alpha_kp1 = intpol3(0, f_0, fp_0, b, f_b, fp_b, lower, upper, PriLev);
      end	% (if)
      %if linalg == 0
      %   alpha_kp1 = intpol2(alpha_k, f_alpha_k, fp_alpha_k,   ...
      %                       alpha_km1, f_alpha_km1, lower, upper, PriLev);
      %elseif LineAlg==1
      %   alpha_kp1 = intpol2(alpha_k, f_alpha_k, fp_alpha_k,   ...
      %                    alpha_km1, f_alpha_km1, lower, upper, PriLev);
      %else
      %   alpha_kp1 = intpol3(alpha_k, f_alpha_k, fp_alpha_k,  alpha_km1, ...
      %                     f_alpha_km1, fp_alpha_km1, lower, upper, PriLev);
      %end
   end	% (if)
   alpha_km1 = alpha_k;
   f_alpha_km1 = f_alpha_k;
   %g_alpha_km1 = g_alpha_k;
   fp_alpha_km1 = fp_alpha_k;
   alpha_k = alpha_kp1;
   k = k + 1;

   if k > MaxIter
      if PriLev > 0
         fprintf('alpha_k: %f, f_alpha_k %f\n',alpha_k,f_alpha_k)
         fprintf('eloLine: Too many PHASE 1 iterations, %d its.\n',k)
      end
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
   fprintf('PHASE 2 === Interval [a, b] = %16.12f %16.12f\n',a,b);
end
fprintf('                     [fa, fb] = %16.12f %16.12f\n',f_a,f_b);

Result2=Result1;
Result1=Result;
Result=[];
pause

while 1
% Extra test 1, end if the interval is too small
%
   %if abs(a-b) <= eps1
   %   if f_a < f_b | (f_a==f_b & a > b)
   %      alpha_k = a;
   %      f_alpha_k = f_a;
   %      g_alpha_k = g_a;
   %      Result=Result1;
   %   else
   %      Result=Result2;
   %      alpha_k = b;
   %      f_alpha_k = f_b;
   %   end
   %   if alpha_k ~= 0
   %      [g_alpha_k, Result] = evalGrad(g, xEval, gType, Result,varargin{:});
   %   end
   %   if PriLev > 0
   %      fprintf('=== End line search, [a b] < eps1 = %16.8e.',eps1);
   %      fprintf(' alpha_k = %16.8e\n',alpha_k);
   %   end
   %   Result.alpha=alpha_k;
   %   Result.alphaVec=alphaVec;
   %   Result.f_alpha=f_alpha_k;
   %   Result.g_alpha=g_alpha_k;
   %   return;
   %end
%
% Compute new alpha using interpolation
%
   if a < b
      % Choose a point in [0,b]
      lower = 0 + tau2*(b-0);
      upper = b - tau3*(b-0);		% Choose alpha in [lower upper]
   else
      % Choose a point in [b,a]
      lower = a + tau2*(b-a);
      upper = b - tau3*(b-a);		% Choose alpha in [lower upper]
   end
   fprintf('PHASE 2 === [lower, upper] = %16.12f %16.12f\n',lower,upper);
   if 1 & NLLS
      % v1 = J_0*p; is computed
      FDIFF = (Result1.r_k-r_0)/a;
      %xprinte(FDIFF,'FDIFF:')
      v2 = (FDIFF-v1)/a;
alpha_k
0.5*(r_0'*r_0)
0.5*(Result1.r_k'*Result1.r_k)
      [alpha_k] = parrm(r_0, v1, v2,a);
      fprintf('parrm-2: %16.14f in%17.14f%17.14f\n',alpha_k,lower,upper);
      alpha_k=min(max(min(lower,upper),alpha_k),max(upper,lower));
   elseif 1 | LineAlg==0
LineAlg
      % Always use f'(0), and two function values, to avoid computing f'(x)
      %alpha_k = intpol2(a, f_a, fp_a, b, f_b, lower, upper, PriLev);
      alpha_k = intpol2(0, f_0, fp_0,   ...
                          a, f_a, lower, upper, PriLev);
      fprintf('intpol2: %16.14f in%17.14f%17.14f\n',alpha_k,lower,upper);
   elseif LineAlg==1
      % Use f'(0) and three function values f'(0),f(a),f(b)
      % MUST REWRITE INTPOL3
      alpha_k = intpol3(a, f_a, fp_a, b, f_b, fp_b, lower, upper, PriLev);
   end	% (if)
   if PriLev > 0
      fprintf('Interpolated alpha, alpha_k = %16.12f\n',alpha_k);
   end
%
% Extra test 2. Check if reduction possible.
%
   %if (a-alpha_k) <= eps2
   %   if PriLev > 0
   %      fprintf('alpha_k=%18.16e; a=%18.16e',alpha_k,a);
   %      fprintf(' END: Too small interval\n');
   %   end
   %   alpha_k = a;
   %   f_alpha_k = f_a;
   %   g_alpha_k = g_a;
   %   if alpha_k ~= 0
   %      [g_alpha_k, Result] = evalGrad(g, xEval, gType, Result,varargin{:});
   %   end
   %   Result.alpha=alpha_k;
   %   Result.alphaVec=alphaVec;
   %   Result.f_alpha=f_alpha_k;
   %   Result.g_alpha=g_alpha_k;
   %   return;
   %end
%
% Compute new f(alpha)
%
   xEval=x+alpha_k*p;
   [f_alpha_k, Result] = evalFunc( f, xEval, fType, Result,varargin{:});
alpha_k
f_alpha_k
pause

   alphaVec=[alphaVec, alpha_k];
   if PriLev > 0
      fprintf('OLD BEST f(a) = %16.8e    a=%16.12f \n',f_a,a);
      fprintf('New  f_alpha_k = %16.8e alpha=%16.12f \n',f_alpha_k,alpha_k);
   end
%
% New alpha better?
%
   if f_alpha_k >= f_a	
alpha_k
a
f_alpha_k
f_a
      % New alpha_k worse, accept the previous point
      alpha = a;		% New alpha_k worse
      fprintf('New alpha_k worse. Stop with %14.12f\n',alpha)
      f_alpha_k = f_a;		%
f_alpha_k
      Result=Result1;
0.5*(Result.r_k'*Result.r_k)
pause
      if alpha_k ~= 0
         [g_alpha_k, Result] = evalGrad(g, xEval, gType, Result,varargin{:});
      end
      Result.alpha=alpha_k;
      Result.alphaVec=alphaVec;
      Result.f_alpha=f_alpha_k;
      Result.g_alpha=g_alpha_k;
      return;
   elseif f_alpha_k >  f_0 + rho*alpha_k*fp_0
      % Above the rho line. Accept point and do one more trial
      disp('Above the rho line. Accept new point')
      b = a;
      f_b = f_b;
      a = alpha_k;              % New alpha the best uptil now
      f_a = f_alpha_k;
      Result2=Result1;
      Result1=Result;
   else				% New alpha_k better than a
      % Now below rho line. Either stop or improve one more step
      % Estimate directed derivative

      fp_alpha_k=(f_alpha_k-f_a)/(alpha_k-a);
      disp('Below rho line. Test slope')
      fp_alpha_k

      %if (alpha_k-a) < 0.3
      %   fp_alpha_k=(f_alpha_k-f_a)/(alpha_k-a);
      %   disp('Below rho line. Test slope')
      %   fp_alpha_k
      %  pause
      %else
      %   fp_alpha_k=Inf;
      %end

      if abs(fp_alpha_k) <= -sigma*fp_0
         if PriLev > -100
            fprintf('f_alpha_k = %16.8e alpha=%16.8e\n',f_alpha_k,alpha_k);
            fprintf('fp_alpha_k = %16.8e fp_0=%16.8e\n',fp_alpha_k,fp_0);
            disp('--- End of line search, F-2 ---');
         end
         if alpha_k ~= 0
            [g_alpha_k, Result] = evalGrad(g, xEval, gType, Result,varargin{:});
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
      % NOT SUFFICIENT.
      if (b-a)*fp_alpha_k >= 0 & ~isinf(fp_alpha_k) % a is new b
         b = a;
         f_b = f_a;
         %g_b = g_a;
         fp_b = fp_a;
         Result2=Result;
      end	% (if)
      a = alpha_k;		% New alpha the best uptil now
      f_a = f_alpha_k;
      fp_a = fp_alpha_k;
      Result1=Result;
      Result=[];
   end	% (else)
end	% (while)

% MODIFICATION LOG:
%
% 980918  hkh  Changed f_min to fLow
% 980919  hkh  Changed to use input structure optParam instead of vector param
% 981005  hkh  Change to more general structure Result as output, including
%              line search on merit functions and penalty functions.
%              Changing function name from linesrch.m to eloLine.m.
%              Changing variable names.
% 981017  hkh  Errors in sending the correct c_k to g (pType==3,4).
%              Changing the logic for the merit function computations
% 981019  hkh  Not setting Result.f_k to f_0 at start, error if alpha=0 and
%              merit function.
% 981027  hkh  Move fLow to optParam.fLow
% 981110  hkh  Change dc_k to cJac in Result and OUT structure
%              Change to alpha in comments and printings
% 001109  hkh  Use LineParam, not optParam, optParam.LineSearch
% 090717  med  Some calculations updated

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
   fprintf('Quad input x0=%14.12f, a=%15.12f, b=%15.12f x1=%15.6e\n',x0,a,b,x1);
   fprintf('Quad func  f0 %14.12f,f1=%15.12f,g0=%15.12f\n',f0,f1,g0);
end
if isinf(f0)
   error('intpol2: Gradient f0 has infinite value')
   return
end
if isinf(g0)
   error('intpol2: Gradient g0 has infinite value')
   return
end
% Transform derivative to [0,1]
g0 = g0 * (x1-x0);
% Let c = (f1 - f0 - g0)
c = f1 - f0 - g0;
% Compute minimum
if c==0
   if PriLev > 3
      fprintf('intpol2: ERROR - computing quadratic interpolation\n')
      fprintf('g0 = %15.7e, c  = %15.7e ',g0,c)
      fprintf('f1 = %15.7e, f0 = %15.7e \n',f1,f0)
   end
   alfa=min(a,b);
   return
else
   z = (-g0/c)/2;
end
% Function value at minium
q_z = f0 + g0*z + c*z^2;
% Compute function value at endpoints in interval [a,b]
% Must transform [a,b] to [0,1] ( = [x0,x1] )
% z_new = (alfa-a)/(b-a) = (?-x0)/x1-x0)
z_a = (a-x0)/(x1-x0);
q_a = f0 + g0*z_a + c*z_a^2;
z_b = (b-x0)/(x1-x0);
q_b = f0 + g0*z_b + c*z_b^2;
if PriLev > 3
   fprintf('Quad minimum  %14.12f, a= %14.12f, b= %14.12f\n',z,z_a,z_b);
   fprintf('Quad func.val %14.12f,fa= %14.12f,fb= %14.12f\n',q_z,q_a,q_b);
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
   fprintf('Cubic input x0=%10.6f, a= %10.6f, b= %10.6f x1=%10.6f\n',x0,a,b,x1);
   fprintf('Cubic func  f0 %10.6f,f1= %10.6f,g0=%10.6f, g1=%10.6f\n',...
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
      fprintf('s, z^3-term is too small; s %12.8e r %12.8e \n',s,r);
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
      fprintf('Cubic minimum 1 %10.6f, z_1= %10.6f\n',c_1,z_1);
      fprintf('Cubic minimum 2 %10.6f, z_2= %10.6f\n',c_2,z_2);
      fprintf('Min in oldCoord %10.6f,    = %10.6f\n',...
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
      fprintf('Interval in z:  a= %10.6f, b= %10.6f\n',z_a,z_b);
      fprintf('Cubic func.val fa= %10.6f,fb= %10.6f\n',c_a,c_b);
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
      fprintf('Best alfa %10.6f\n',alfa);
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
function [f_k, Result] = evalFunc( f, x, fType, Result, varargin)
% ============================================================

nargin;
switch fType
case 1
   f_k = nlp_f(x,varargin{:});
case 2
   r_k = nlp_r(x , varargin{:} );     % Function value
   f_k = 0.5 * (r_k' * r_k);
   Result.r_k = r_k;
case 3
   [f_k Result] = con_fm(x, Result, varargin{:});
   % Set in con_fm: Result.f_k Result.c_k Result.Ax
case 4
   f_k = feval(f,x,varargin{:});
case {5,6}
   r_k = feval(f,x,varargin{:});
   f_k = 0.5 * (r_k' * r_k);
   Result.r_k = r_k;
case {7,8}
   [f_k Result] = feval(f, x, Result, varargin{:});
   % Set in con_fm: Result.f_k Result.c_k Result.Ax
end

% ============================================================
function [g_k, Result] = evalGrad( g, x, gType, Result, varargin)
% ============================================================

nargin;
switch gType
case 1
   g_k = nlp_g(x, varargin{:});
case 2
   J_k = nlp_J(x , varargin{:} );
   g_k = J_k' * Result.r_k;
   Result.J_k = J_k;
case 3
   [g_k Result] = con_gm(x, Result, varargin{:});
   % Set in con_gm: Result.g_k Result.cJac
case 4
   g_k = feval(g, x, varargin{:});
case {5,6}
   J_k = feval(g, x, varargin{:});
   g_k = J_k' * Result.r_k;
   Result.J_k = J_k;
case {7,8}
   [g_k Result] = feval(g, x, Result, varargin{:});
   % Set in con_gm: Result.g_k Result.cJac
end

%---------------------------------

function[alpha] = parrm(v0, v1, v2, alpha_in)

v1n = norm(v1)^2;
v2n = norm(v2)^2;

if v2n==0
   alpha=1000;
   return
end

a1 = (1.5*(v1'*v2))/v2n;
a2 = (v0'*v2+.5*v1n)/v2n;
a3 = (0.5*(v0'*v1))/v2n;

aa(1) = 1;
aa(2) = a1;
aa(3) = a2;
aa(4) = a3;
rr = roots(aa);
rr = rr(find(imag(rr)== 0));
rr = rr(find((rr) > 0));

alpha=min(rr);

% Now we have the intersting root, rr.

function [zz]= Newpnt(zzz)

%This function halves the steplength alfin until the actual
%decrease of the objective function is at least tau*the predicted.

n=length(x);
v0 = f;
% J is not defined
%v1 = J*p;
F0 = f;

tau=0.1;
xny=x+p;

% Compute f

F2 = [fny;mu*xny];
alf2= alfa;
sec = 0;

while (Fsq0-Fsq)<tau*(Fsq0-0.5*norm([f0;mu*x]+alfa*[Jac0;mu*eye(n)]*p)^2)
  if sec == 0 % Interpolate one gradient and two function values
	FDIFF = F2-F0;
	v2 = FDIFF-v1;
	[alfa] = parrm(v0, v1, v2,alfa);
  	xny=x+alfa*p;
 %Compute f
  	i=i+1;
	alf1 = alfa;
  end
  if sec == 1
	U1 = (F1-F0) / alf1;
	U2 = ((F2-F1)/(alf2-alf1) - (F1-F0)/alf1) / alf2;
	v1 = U1-alf1*U2;
	v2 = U2;
 %Compute f
	alf2 = alf1;
	F2 = F1;
	alf1 = alfa;
	i = i + 1;
  end
  sec = 1;
end
eval=i;
stplen=alfa;
