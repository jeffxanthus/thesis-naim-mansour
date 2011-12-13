% glcSolve.m
% 
% Solves general constrained mixed integer global optimization problems.
%
% glcSolve.m implements the algorithm DIRECT by Donald R. Jones presented
% in the paper "DIRECT", Encyclopedia of Optimization, Kluwer Academic
% Publishers, 2001.
%
% glcSolve solves problems of the form:
%
% min   f(x)
%  x
% s/t   x_L <=   x  <= x_U
%       b_L <= A x  <= b_U
%       c_L <= c(x) <= c_U
%       x(i) integer, for i in I
%
% Recommendation: Put the integers as the first variables !
% Put low range integers before large range integers
% Linear constraints are specially treated
% Equality constraints are added as penalties to the objective.
% Weights are computed automatically, assuimg f(x) scaled to be roughly 1
% at optimum. Otherwise scale f(x)
%
% See below "USAGE" on how to create the Prob structure and do the call
%
% Calling syntax:
%
% function Result = glcSolve(Prob, varargin )
%
% INPUT PARAMETERS
%
% Prob    Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   FUNCS.f   The routine to compute the function, given as a string, say GLCF
%   FUNCS.c   The routine to compute the nonlinear constraint, say GLCC
%             A call to tomFiles could be used to set the names into
%             the Prob struct:
%             Prob = tomFiles(Prob,'GLCF',[],[],'GLCC');
%   x_L       Lower bounds for each element in x.
%   x_U       Upper bounds for each element in x.
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
%   PriLevOpt Print Level 
%             0 = silent. 1 = some printing. 2 = print each iteration
%   WarmStart If true, >0, glcSolve reads the output from the last run
%             from the mat-file glcSave.mat, and continues from the last run.
%
% optParam    Structure in Prob, Prob.optParam. 
%             Defines optimization parameters. Fields used:
%   IterPrint Print one line each iteration
%   cTol      Nonlinear constraint tolerance
%   MaxIter   Maximal number of iterations, default 50.
%   MaxFunc   Maximal number of function evaluations, default 200 (roughly).
%   EpsGlob   Global/local weight parameter, default 1E-4.
%   fGoal     Goal for function value, if empty not used
%   eps_f     Relative accuracy for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol , if fGoal \=0
%             Stop if abs(f-fGoal) <= fTol , if fGoal ==0
%   eps_x     Convergence tolerance in x. All possible rectangles are 
%             less than this tolerance (scaled to (0,1) )
%             See the output field maxTri.
%
% MIP         Structure in Prob, Prob.MIP.
%             Defines integer optimization parameters. Fields used:
%   IntVars   Set of integer variables, default Integers=[].
% GO          Structure in Prob, Prob.GO.
%   fEqual    All points with function values within tolerance fEqual are 
%             considered to be global minima and returned
%   LinWeight RateOfChange = LinWeight*|a(i,:)| for linear constraints.
%             Balance between linear and nonlinear constraints.
%             Default value 0.1.
%       
% OUTPUT PARAMETERS
%
% Result    Structure with results from optimization
%  x_k      Matrix with optimal points as columns.
%  f_k      The best function value found so far
%  c_k      Nonlinear constraints values at x_k
%  Iter     Number of iterations
%  FuncEv   Number of function evaluations
%  ExitText Text string giving ExitFlag and Inform information
%
% To make a warm start possible, glcSolve saves the following information in
% the file glcSave.mat:
%   Name      Name of the problem. Used for security if doing warm start
%   C         Matrix with all rectangle centerpoints, in [0,1]-space.
%   F         Vector with function values.
%   T         T(i) is the number of times rectangle i has been trisected.
%   D         Vector with distances from centerpoint to the vertices.
%   G         Matrix with constraint values for each point.
%   I_L       I_L(i,j) is the lower bound for rect. j in integer dim. I(i)
%   I_U       I_U(i,j) is the upper bound for rect. j in integer dim. I(i)
%   s_0       s_0 is used as s(0)
%   s         s(j) is the sum of observed rates of change for constraint j.
%   t         t(i) is the total # splits along dimension i.
%   ignoreidx Rectangles to be ignored in the rect. selection proceedure.
%   feasible  Flag indicating if a feasible point has been found.
%   Split     Split(i,j) = # splits along dimension i of rectangle j
%   glcfMin   Best function value found at a feasible point.
%   fMinIdx   Indices of the currently best points
%   fMinEQ     sum(abs(infeasibilities)) for minimum points, 0 if no equalities
%
% USAGE:
%
% The function GLCF is best written as
%     function f = GLCF(x, Prob) 
% Then any information, say u and W is easily sent to GLCF (and GLCC) using 
% the Prob structure. See the example below. 
%
% Assume bounds x_L and x_U are initialized, as well as the linear constraint
% matrix A, lower and upper bounds, b_L and b_U on the linear constraints,
% and lower and upper bounds c_L and c_U on the nonlinear constraints
% (Put [] if all bounds are inf or -inf). Use the TOMLAB Quick format:
%
%      The name of the problem is set to "GLCF Test"
%
%      Prob   = glcAssign('GLCF',x_L,x_U,'GLCF Test',A,b_L,b_U,'GLCC',c_L,c_U);
%      Prob.user.u = u; Prob.user.W=W;    % example of extra user data
%
%      % Default values are now set for PriLevOpt, and structure optParam
%      % To change a value, an example is shown on the next line
%      Prob.optParam.MaxFunc = 500; % Change max number of function evaluations 
%
%      If there are integer variables, they may be set as additional input
%      to glcAssign, or directly as the next line shows:
%      Prob.MIP.IntVars = [1 3];  % 1st and third variables are integers
%
%      Result = glcSolve(Prob);
%      PrintResult(Result);
%           
% The user function GLCF is written as
%
%      function f = GLCF(x, Prob) 
%      u = Prob.user.u; W = Prob.user.W;
%      f = "some function of x, u and W"
%
% It is also possible to use the function format
%      function f = GLCF(x) 
% or sending the extra parameters as additional input 
%
%      Result = glcSolve(Prob,u,W);
%
%      function f = GLCF(x,Prob,u,W) 
%
% NOTE! If additional parameters are sent, Prob must be the second input 
% parameter to GLCF (and GLCC)
%
% The user function GLCC, computing the nonlinear constraints, is written as
%
%      function c = GLCC(x, Prob) 
%      u = Prob.user.u; W = Prob.user.W;
%      c = "some vector function of x, V and W"
%
% To make a restart, just set the restart flag, and call glcSolve once again:
%
%      Prob.WarmStart = 1;
%      Result = glcSolve(Prob);   % Assuming no extra parameters in the call
%      PrintResult(Result);
%
% If the maximum number of iterations have been reached, it must be increased:
%      Prob.optParam.MaxIter = 1000;
%
% To change the number of function evaluations for each call to glcFast:
%      Prob.optParam.MaxFunc = 500;

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomopt.com
% Copyright (c) 1999-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Feb 15, 1999.   Last modified Aug 14, 2006.

function Result = glcSolve(Prob, varargin )

if nargin < 1
   error('glcSolve needs input structure Prob');
end

solvType=checkType('glc');

Prob=ProbCheck(Prob,'glcSolve',solvType);

if isempty(Prob.x_L) | isempty(Prob.x_U)
   disp('glcSolve requires both lower and upper variable bounds');
   Result.ExitFlag = 1;
   Result.ExitText = 'glcSolve requires both lower and upper variable bounds';
   Result=endSolve(Prob,Result);
   return;
end

Prob=iniSolve(Prob,checkType('glc'),0,0);

DEBUG=0;

PriLev    = Prob.PriLevOpt;          % Print level
MaxFunc   = Prob.optParam.MaxFunc;   % Maximum number of function evaluations
MaxIter   = Prob.optParam.MaxIter;   % Number of iterations
EpsGlob   = Prob.optParam.EpsGlob;   % global/local weight parameter. 
cTol      = Prob.optParam.cTol;      % Constraint feasibility tolerance
bTol      = Prob.optParam.bTol;      % Linear Constraint feasibility tolerance
IterPrint = Prob.optParam.IterPrint; % Print short information each iteration
fGoal     = Prob.optParam.fGoal;     % Goal for f(x).
fTol      = Prob.optParam.eps_f;     % Relative tolerance for fGoal
xTol      = Prob.optParam.eps_x;     % Tolerance for rectangle sizes
                                     % (scaled to (0,1) )
%MaxIter = 5000
%MaxFunc = 5000
%%MaxFunc = 15
%IterPrint = 1


glcfOld = Inf;

if isfield(Prob.MIP,'IntVars')
   IntVars   = Prob.MIP.IntVars(:);     % Set of integer variables.
else
   IntVars = [];
end
if isfield(Prob.GO,'LinWeight')
   LinWeight   = Prob.GO.LinWeight;     % Set of integer variables.
else
   LinWeight = [];
end
if isempty(LinWeight), LinWeight = 0.1; end

if isfield(Prob.GO,'fEqual')
   fEqual   = Prob.GO.fEqual;     % Tolerance for points being equal
else
   fEqual = [];
end
if isempty(fEqual), fEqual = 1E-10; end

betaFac  = 10;
nFunc    = 0;
convflag = 0;
AvIter   = 100;

Result                 = ResultDef(Prob);
Result.Solver          = 'glcSolve';
Result.SolverAlgorithm = 'DIRECT - Lipschitzian Optimization';


% Pick up input parameters from the Prob structure:
x_L = Prob.x_L(:);   % Lower bounds
x_U = Prob.x_U(:);   % Upper bounds
A   = Prob.A;        % Linear constraint matrix
b_L = Prob.b_L(:);   % Lower bounds, linear constraints
b_U = Prob.b_U(:);   % Upper bounds, linear constraints
c_L = Prob.c_L(:);   % Lower bounds, nonlinear constraints
c_U = Prob.c_U(:);   % Upper bounds, nonlinear constraints

n     = length(x_L);  % Problem dimension
nFunc = 0;            % Function evaluation counter

isMIP  = ~isempty(IntVars); % isMIP is true if there are integer variables
I_logic          = zeros(n,1); % Logic set for integer/continuous variables
I_logic(IntVars) = 1;

% Index sets for linear and nonlinear constraints
beta = [];
if isempty(b_L)
   leq     = [];
   b_L_idx = [];
   if isempty(b_U)
      b_U_idx = [];
   else
      b_U_idx = find(isfinite(b_U));
   end
elseif isempty(b_U)
   leq     = [];
   b_U_idx = [];
   b_L_idx = find(isfinite(b_L));
else
   leq     = b_L == b_U;
   b_L_idx = find(isfinite(b_L) & ~leq);
   b_U_idx = find(isfinite(b_U) & ~leq);
   leq     = find(leq);
   if ~isempty(leq)
      v       = b_L(leq);
      v(v==0) = 1;
      % Weight on infeasibility addition
      beta    = betaFac./v;
   end
end
if isempty(c_L)
   nleq    = [];
   c_L_idx = [];
   if isempty(c_U)
      c_U_idx = [];
   else
      c_U_idx = find(isfinite(c_U));
   end
elseif isempty(c_U)
   nleq    = [];
   c_U_idx = [];
   c_L_idx = find(isfinite(c_L));
else
   nleq    = c_L == c_U;
   c_L_idx = find(isfinite(c_L) & ~nleq);
   c_U_idx = find(isfinite(c_U) & ~nleq);
   nleq    = find(nleq);
   if ~isempty(nleq)
      v       = c_L(nleq);
      v(v==0) = 1;
      % Weight on infeasibility addition
      beta    = [beta;betaFac./v];
   end
end


%b_L_idx = find(~isempty(b_L) & isfinite(b_L));
%b_U_idx = find(~isempty(b_U) & isfinite(b_U));
%c_L_idx = find(~isempty(c_L) & isfinite(c_L));
%c_U_idx = find(~isempty(c_U) & isfinite(c_U));

% Check if there is a chance to find a feasible point
if any(x_L > x_U)
   fprintf('\n\n');
   fprintf('Error in glcSolve, upper variable bounds below lower bounds\n');
   Result.ExitFlag = 2;
   Result.ExitText =str2mat('glcSolve solves box-bounded problems' ...
            ,'Found some upper bound below lower bounds');
   Result=endSolve(Prob,Result);
   return;
end
x_D   = x_U-x_L;
nFix  = sum(x_D == 0);
nCont = sum(x_D > 0 & ~I_logic);

%
%  STEP 1, INITIALIZATION
%

if Prob.WarmStart
   % Restart with values from previous run.

   load('glcSave.mat','Name','C','F','T','D','G','I_L','I_U','s_0','s','t',...
   'ignoreidx','feasible','Split','glcfMin','fMinIdx','fMinEQ');
   Name1 = Prob.Name;               % Name for the problem

   if strcmp(Name1,Name)
      LinI  = length(b_L_idx) + length(b_U_idx);% # of linear constraints
      NLinI = length(c_L_idx) + length(c_U_idx);% # of nonlinear constraints
      mI    = LinI + NLinI;                     % Number of constraints
      LinE  = length(leq);
      NLinE = length(nleq);
      mE    = LinE + NLinE;
      m     = mI + mE;
      if m > 0
         % Place all constraints in one vector gx with upper bounds g_U
         g_U = [-b_L(b_L_idx);b_U(b_U_idx);-c_L(c_L_idx);c_U(c_U_idx); ...
               b_L(leq);c_L(nleq)];
         g_T = [bTol*max(1,abs(b_L(b_L_idx)));bTol*max(1,abs(b_U(b_U_idx)));...
                cTol*max(1,abs(c_L(c_L_idx)));cTol*max(1,abs(c_U(c_U_idx)));...
                bTol*max(1,abs(b_L(leq)));cTol*max(1,abs(c_L(nleq)))];
         G   = [G, zeros(m,MaxFunc + 20)];
      else
         rE = 0;
      end
      AvIter = AvIter-min(AvIter,nFunc); % Either 0, or iterations left
      AvFunc = 0;
      nFuncOld = size(F,2);
      % T(i), number of times rectangle i has been trisected
      T = [T;zeros(MaxFunc-length(T),1)]; 

      if PriLev > 0
         fprintf('\n');
         fprintf('Restarting with %d sampled points from previous run\n', ...
         size(F,2));
      end
      if ~isinf(glcfMin) 
         % All points i with F(i)=glcfMin, and feasible
         xBest = tomsol(9,x_L,C(:,fMinIdx),x_D);    
      else
         xBest = [];
      end
   else
      Prob.WarmStart = 0;
      if PriLev >= -1000
         fprintf('Previous run was with Problem %s\n',Name);
         fprintf('This run is with Problem %s\n',Name1);
         fprintf('Impossible to do restart.\n');
         fprintf('Maybe there exists several files glcSave.mat?\n');
      end
   end
end % if WarmStart

if ~Prob.WarmStart
   Name = deblank(Prob.Name);  % Problem name
   % No restart, set first point to center of the unit hypercube.
   nFuncOld = 0;

   % SAMPLE THE CENTERPOINT OF THE ENTIRE SPACE.
   C = ones(n,1)./2;    % Matrix with all rectangle centerpoints.
   % All C_coordinates refers to the n-dimensional hypercube. 
   Split = zeros(n,1); % Split(i,j) is the number of times rectangle j has 
                       % been split along dimension i. Split(i,j) is set to
                       % Inf if rect j has sidelength 0 for integer variable i. 

   % IF ALL VARIABLES ARE INTEGERS, THEN IT IS POSSIBLE FOR A RECTANGLE
   % TO BE REDUCED TO A SINGLE POINT. IF THIS HAPPENS, THE RECTANGLE SHOULD
   % BE IGNORED IN THE RECTANGLE SELECTION PROCEEDURE.
   ignoreidx = []; 

   % Fixed variables should not be considered for the splitting procedure
   idx = find(x_L==x_U);
   % Eliminates risk choosing these dims as splitting dims for any rectangle 
   if length(idx) > 0
      Split(idx,1) = Inf; 
      if length(idx)==n
         ignoreidx = [ignoreidx 1];
      end
   end

   if isMIP
      % I_L(i,j) is the lower bound for rectangle j in integer dimension I(i).
      I_L =  ceil(x_L(IntVars)); 
      % I_U(i,j) is the upper bound for rectangle j in integer dimension I(i).
      I_U = floor(x_U(IntVars)); 
      IL = (I_U-I_L);
      if any(IL < 0)
         disp('Error in glcSolve, empty domain for integer variables:')
         tmpidx = find(IL<0);
         IL(tmpidx)
         Result.ExitFlag = 8;
         Result.ExitText = 'Empty domain for integer variables';
         Result=endSolve(Prob,Result);
         return;
      end
      x_mid  = floor( (I_U+I_L)/2 );
      %x_mid  = floor( (x_U(I)+x_L(I))/2 );
      C(IntVars,1) = (x_mid-x_L(IntVars))./(max(x_D(IntVars),1E-20));
      tmpidx    = find(IL(:,1)==0);
      % Eliminate risk choosing these dims as splitting dim for any rectangle
      if ~isempty(tmpidx)
         Split(IntVars(tmpidx),1) = Inf; 
         if length(tmpidx)==n
            ignoreidx = [ignoreidx 1];
         end 
      end
   else
      I_L = [];
      I_U = [];
   end

   % Transform C to original search space
   %x     = x_L + C.*x_D;  
   x    = tomsol(9,x_L, C,x_D); 
   f     = nlp_f(x, Prob, varargin{:});  % Function value at x
   nFunc = nFunc+1;                      % Number of function evaluations
   T = zeros(MaxFunc,1); % T(i), number of times rectangle i has been trisected
   D = sqrt(n-nFix)/2; % Vector with distances from centerpoint to the vertices

   % IF THE CENTER IS FEASIBLE, SET x_min EQUAL TO THE CENTERPOINT AND
   % glcfMin EQUAL TO THE OBJECTIVE FUNCTION VALUE AT THIS POINT.

   LinI  = length(b_L_idx) + length(b_U_idx);% # of linear constraints
   NLinI = length(c_L_idx) + length(c_U_idx);% # of nonlinear constraints
   mI    = LinI + NLinI;                     % Number of constraints
   LinE  = length(leq);
   NLinE = length(nleq);
   mE    = LinE + NLinE;
   m     = mI + mE;
   if LinI+LinE > 0
      Ax = A*x;
   else
      Ax=[];
   end   
   if NLinI+NLinE > 0 
      cx = nlp_c(x, Prob, varargin{:});
      cx = cx(:);
   else
      cx=[];
   end
   %if LinE > 0
   %   % Use penalty approach for linear equality constraints
   %   %f = f + sum(abs(Ax(leq,:)-b_L(leq)))/bTol;
   %   %f = f + beta*sum(abs(Ax(leq,:)-b_L(leq)));
   %end
   %if NLinE > 0
   %   % Use penalty approach for nonlinear equality constraints
   %   %f = f + sum(abs(cx(nleq,:)-c_L(nleq)))/cTol;
   %   %f = f + beta*sum(abs(cx(nleq,:)-c_L(nleq)));
   %end
   F = f;  % Vector with function values

   if m > 0
      % Place all constraints in one vector gx with upper bounds g_U
      g_U = [-b_L(b_L_idx);b_U(b_U_idx);-c_L(c_L_idx);c_U(c_U_idx); ...
              b_L(leq);c_L(nleq)];
      g_T = [bTol*max(1,abs(b_L(b_L_idx)));bTol*max(1,abs(b_U(b_U_idx)));...
             cTol*max(1,abs(c_L(c_L_idx)));cTol*max(1,abs(c_U(c_U_idx)));...
             bTol*max(1,abs(b_L(leq)));cTol*max(1,abs(c_L(nleq)))];
      gx  = [-Ax(b_L_idx);Ax(b_U_idx);-cx(c_L_idx);cx(c_U_idx); ...
              Ax(leq);cx(nleq)]; 
      tmp = gx - g_U;
      if all( tmp(1:mI) <= g_T(1:mI) ) % Initial point is feasible ?
         feasible = 1; % feasible point has been found
      elseif all( tmp(1:LinI) <= g_T(1:LinI) ) % Linear constraints feasible?
         feasible = 0;
      else
         feasible = -1;
      end
      G      = zeros(m,MaxFunc + 20);
      G(:,1) = tmp; % Subtract g_U to get g(x)<=0 (or g(x) = 0 for ineq.)
      % G is a matrix with constraint values for each point.
      % G(i,j) is the value of constraint i at point j.
   else
      G = [];
      feasible = 1;
   end
   fMinIdx = 1;
   xBest   = x;
   if mE > 0
      fMinEQ = sum(beta.*abs(G(mI+1:m,nFunc+nFuncOld)));
   else
      fMinEQ = 0;
   end
   if feasible > 0
      glcfMin = f+fMinEQ;
   else
      glcfMin = NaN;
   end

   % SET s(j)=0 FOR j=1,2,...,m.
   if m > 0
      s_0 = 0;          % Used as s(0).
      % s(j) is the sum of observed rates of change for constraint j.
      s   = zeros(m,1); 
      if LinI + LinE > 0
         z = sqrt(sum(A.^2'))';
      end
      if LinI > 0
         s(1:LinI) = LinWeight*[z(b_L_idx);z(b_U_idx)];
      end
      if LinE > 0
         s(mI+1:mI+LinE) = LinWeight*z(leq);
      end
   else
      s_0 = [];
      s = [];
      rE = 0;
   end

   AvFunc = 0;

   % SET t(i)=0 FOR i=1,2,...,n.
   % t(i) is the number of times a rectangle has been split along dimension i.
   t = zeros(n,1); 
end % if not restart

%if m == 0
%   cGD = zeros(1,MaxFunc); 
%end

Iter = 0; % Iteration counter
convflag  = 0;

while Iter < MaxIter & nFunc < MaxFunc & convflag == 0 

   Iter = Iter+1;

   %
   %  STEP 2, SELECT RECTANGLES 
   %

   if length(ignoreidx)==length(F) % If all rectangles are fathomed
      break
   end

   % COMPUTE c(j) VALUES USING THE CURRENT VALUES OF s_0 AND s(j), j=1,2,...,m  
   %if m > 0
   %   c = s_0./(max(s,1E-5));
   %   cGD = (c(1:mI)'*max(G(1:mI,1:nFunc+nFuncOld)-cTol,0))./D;
   %   %cGD = max(c*ones(1,nFunc+nFuncOld).*max(G(:,1:nFunc+nFuncOld)-...
   %   %          cTol,0))./D;
   %end

   S = []; % The set of rectangles selected for trisection

   ixD = find(D == 0);
   D(ixD) = 1;
   if feasible < 0
      % Not feasible with respect to the linear constraints
      %if m > 0
         c = s_0./(max(s,1E-5));
         cGD = (c(1:LinI)'*max(G(1:LinI,1:nFunc+nFuncOld)-bTol,0))./D;
         %cGD = max(c*ones(1,nFunc+nFuncOld).*max(G(:,1:nFunc+nFuncOld)-...
         %          cTol,0))./D;
      %end
      % IF A FEASIBLE TRIANGLE HAS NOT BEEN FOUND, SELECT THE RECTANGLE THAT
      % MINIMIZES THE RATE OF CHANGE REQUIRED TO BRING THE WEIGHTED CONSTRAINT
      % VIOLATIONS TO ZERO.
      % How to treat the integer variables????
      cGD(ixD) = Inf;

      if ~isempty(ignoreidx)
         ix0   =ones(length(cGD),1);
         ix0(ignoreidx)=0;
         ix    = find(ix0);
         [a b] = min( cGD(ix) );
         if ~isfinite(a)
            disp(' No feasible integer point exist')
            S = [];
            break;
         end
         S = ix(b);
      else
         [a S] = min( cGD );
      end
   elseif feasible < 1
      %if m > 0
         c = s_0./(max(s,1E-5));
         cGD = (c(1:mI)'*max(G(1:mI,1:nFunc+nFuncOld)-cTol,0))./D;
         %cGD = max(c*ones(1,nFunc+nFuncOld).*max(G(:,1:nFunc+nFuncOld)-...
         %          cTol,0))./D;
      %end
      % IF A FEASIBLE TRIANGLE HAS NOT BEEN FOUND, SELECT THE RECTANGLE THAT
      % MINIMIZES THE RATE OF CHANGE REQUIRED TO BRING THE WEIGHTED CONSTRAINT
      % VIOLATIONS TO ZERO.
      % How to treat the integer variables????
      cGD(ixD) = Inf;

      if ~isempty(ignoreidx)
         ix0   =ones(length(cGD),1);
         ix0(ignoreidx)=0;
         ix    = find(ix0);
         [a b] = min( cGD(ix) );
         if ~isfinite(a)
            disp(' No feasible integer point exist')
            S = [];
            break;
         end
         S = ix(b);
      else
         [a S] = min( cGD );
      end
      if m > 10000
      %if NLinI > 1
         cGD = max(c*ones(1,nFunc+nFuncOld).*max(G(:,1:nFunc+nFuncOld)-...
                   cTol,0))./D;
         if ~isempty(ignoreidx)
            [a b] = min( cGD(ix) );
            b = ix(b);
            if S~=b
               S = [S,b];
               S
            end
         else
            [a b] = min( cGD );
            if S~=b
               S = [S,b];
               S
            end
         end
      end
      if NLinI > 10000  % Do not use this idea
         % Add extra rectangles when having several nonlinear constraints
         for i = 1:NLinI
             j = LinI + i;
             cGD = max(G(j,1:nFunc+nFuncOld),0)./D;
             iy = cGD > 0;
             if ~isempty(ignoreidx)
                iy(ignoreidx) = 0;
             end
             iy = find(iy);
             if ~isempty(iy)
                [a b] = min( cGD(iy) );
                b = iy(b);
                if ~any(b==S)
                   S = [S,b];
                else
                   b = find( abs(cGD(iy) - a) < 1E-12 );
                   if length(b) == 1, break; end
                   b = iy(b);
                   for k = 1:length(b)
                       if ~any(b(k)==S)
                          S = [S,b(k)];
                          break;
                       end
                   end
                end
             end
         end
         S
      end

      % [a b] = min( cGD );
      % % How to treat the integer variables????
      % if ~isempty(ignoreidx)
      %    cGDtmp = cGD;
      %    % rectangle b is fathomed i.e. b is a single point and can't be split
      %    while ismember(b,ignoreidx) 
      %       cGDtmp(b) = Inf;
      %       [a b] = min( cGDtmp );
      %       if ~isfinite(a)
      %          disp(' No feasible integer point exist')
      %          S = [];
      %          break;
      %       end
      %    end
      %    S = b;
      % else
      %    S = b;
      % end   
   else
      % ON THE OTHER HAND, IF A FEASIBLE POINT HAS BEEN FOUND, IDENTIFY THE
      % SET OF RECTANGLES THAT PARTICIPATE IN THE LOWER ENVELOPE. LET S BE
      % THE SET OF SELECTED RECTANGLES.
      if m > 0
         c = s_0./(max(s,1E-5));
         % Feasible w.r.t. inequalties - Include equalities into cGD
         cGD = zeros(nFunc+nFuncOld,1);
    
         if mI == 0
            cGD = (c'*abs(G(:,1:nFunc+nFuncOld)))./D; 
         elseif mE == 0
            cGD = (c'*max(G(:,1:nFunc+nFuncOld)-cTol,0))./D; 
         else
            cGD = (c(1:mI)'*max(G(1:mI,1:nFunc+nFuncOld)-cTol,0))./D + ...
                  (c(mI+1:m)'*abs(G(mI+1:m,1:nFunc+nFuncOld)))./D; 
         end
         %cGD = max(c*ones(1,nFunc+nFuncOld).*max(G(:,1:nFunc+nFuncOld)-...
         %          cTol,0))./D;
      end
      Epsilon = max(EpsGlob*abs(glcfMin),1E-8);

      f_star = glcfMin-Epsilon;
      %r_previous=Inf;
      r_previous=[];
      idx1=[];
      while 1
         f_star_max = -Inf;
         if m > 0
            h      = max((F-f_star),0)./D + cGD;
         else
            h      = max((F-f_star),0)./D;
         end
         h(ixD) = Inf;
%if nFunc > 247 & Iter > 21
%   [f_star F(247) D(247) h(247)]
%   pause
%end

         tmpset1=ones(1,length(F));
         tmpset1(idx1)      = 0; % not consider rects put in S in previous turn
         tmpset1(ignoreidx) = 0; % not consider rects being fathomed
         tmpidx=find(tmpset1);
 
         if ~isempty(tmpidx)
            h_min  = min(h(tmpidx));
            idx1   = tmpidx(find(h(tmpidx)==h_min));
         else
            idx1       = [];
         end

 
         if length(idx1) > 1
            if any(f_star > F(idx1))
               % Choose rectangle being constant for the lowest value of f_star
               [a b] = max(f_star-F(idx1));
               r = idx1(b);
            else
               % Choose the rectangle with the flatest slope.
               [a b] = max(D(idx1));
               r = idx1(b);
            end
         else
            r = idx1;
         end
         S = [S setdiff(idx1,S)];
%xprinti(S,'S:',4,23);

         if f_star > F(r) % We must move horisontally to the left
            f_star = F(r);
         end

         % Compute the intersection of the line y(x) = slope*x + const,
         % with all curves not in [idx;r_previous]
         slope = -1/D(r);
         const = h(r)-slope*f_star;
 
         % Try speed up search for f_star_max just checking relevant rectangles

         %idxspeed = find( (cGD<cGD(r)) | (D>D(r)) );
         %idxspeed = setdiff(idxspeed,[idx1 r_previous]);  

         if m > 0
            ix0 =  (cGD<cGD(r)) | (D>D(r)) ;
         else
            ix0 =  D>D(r);
         end
         ix0(r_previous) = 0;
         ix0(idx1)       = 0;
         ix = find(ix0);

         %if isempty(ix) & isempty(idxspeed)
         %else
         %if ~all(ix == idxspeed)
         %   keyboard
         %end
         %end
         %f_star_m = f_star_max;

         if ~isempty(ix)
            if m > 0
               tmp1 = (cGD(ix)-const)/slope;
            else
               tmp1 = (-const/slope)*ones(1,length(ix));
            end

            iy = find(tmp1 < f_star & tmp1 >= F(ix));
            f_star_max = max([f_star_max,tmp1(iy)]);
         end

         sl_hat = -1./D(ix);
         iz = find(sl_hat > slope);
         if m > 0
            con_hat = cGD(ix(iz)) - sl_hat(iz).*F(ix(iz));
         else
            con_hat = -sl_hat(iz).*F(ix(iz));
         end
         tmp2 = (con_hat-const)./(slope-sl_hat(iz));
         iv = find(tmp2 < f_star);
         if ~isempty(iv)
            f_star_max = max([f_star_max,tmp2(iv)]);
         end

         %for dummy=1:length(idxspeed)
         %    i = idxspeed(dummy);
         %    % First assume that h(i) is constant i.e. max(F(i)-tmp1,0)=0
         %    tmp1 = (cGD(i)-const)/slope; % same as tmp1=-D(r)*(cGD(i)-const);
         %    if (tmp1 < f_star) & (tmp1 >= F(i) )
         %       % intersection between y(x) and h(i) where h(i) is constant.
         %       if tmp1 > f_star_max
         %          f_star_max = tmp1;
         %       end
         %    end
         %    %else % HKH-Maybe this should be an end statement !!!!!!!!
         %       % Assume h(i) is not constant i.e. max(F(i)-tmp2,0)=F(i)-tmp2
         %       slope_hat = -1/D(i);
         %       % Else, intersection not occurs or occur for values >= f_star
         %       if slope_hat > slope 
         %          %const_hat = (F(i)-tmp1)/D(i) + cGD(i) - slope_hat*f_star;
         %          const_hat = cGD(i) - slope_hat*F(i);
         %          tmp2 = (const_hat-const)/(slope-slope_hat);
         %          if tmp2 < f_star & tmp2 > f_star_max
         %             f_star_max = tmp2;
         %          end
         %       end
         %    %end
         %end % dummy=1:length(idxspeed)   

         %if f_star_max ~= f_star_m
         %   f_star_max
         %   f_star_m
         %   keyboard
         %end
 
         if isfinite(f_star_max) 
            f_star = f_star_max;
            r_previous = r;
         else
            % if curve r never intersected with another one
            break;
         end
         % DO NOT DOUBLE-COUNT ANY RECTANGLE !!!!!!   
      end % while 1  
   end % if ~feasible
   D(ixD) = 0; % Reset lengths to 0


   % -----------------------------------
   % STEP 3, CHOOSE ANY RECTANGLE IN S. 
   % -----------------------------------


   %if DEBUG      
   %   DDDD=find(all(Split==Inf));
   %   Disaster=intersect(S,DDDD);
   %   if length(Disaster>0)
   %      Disaster
   %      S
   %      ignoreidx
   %      pause
   %   end
   %end

%[Iter,length(S)]
%xprinti(S,'S:',4,23);

   for dummy=1:length(S) % for each rectangle in S
       r = S(dummy);

if DEBUG
r
end

       % ---------------------------------------
       % STEP 4, TRISECT AND SAMPLE RECTANGLE r
       % ---------------------------------------

       % CHOOSE A SPLITTING DIMENSION BY IDENTIFYING THE SET OF LONG SIDES OF
       % RECTANGLE r AND THEN CHOOSE THE LONG SIDE WITH THE SMALLEST t(i) VALUE.
       % IF MORE THAN ONE SIDE IS TIED FOR THE SMALLEST t(i) VALUE, CHOOSE THE 
       % ONE WITH THE LOWEST DIMENSIONAL INDEX.
       idx2 = find(Split(:,r)==min(Split(:,r)));
       if length(idx2) > 1
          idx3 = find(t(idx2)==min(t(idx2)));
          i = idx2(idx3(1));
          %i = idx2(idx3(end));
       else
          i = idx2;
       end
       %if Split(i,r) == 0 & nFunc > 1000
       %   keyboard
       %end

       %if DEBUG %I_L(i,r)==I_U(i,r)
       %   fprintf('\n\n Error in glcSolve, I_L(i,r)==I_U(i,r) !!!!!! \n');
       %   r
       %   i
       %   Split(:,r)
       %   C(:,r)
       %   pause
       %end


       % Updates
       t(i) = t(i)+1;

       %T(r) = T(r)+1;

       %% Update D(r)
       %j = mod(T(r),n-nFix);
       %k = (T(r)-j)/(n-nFix);
       %D(r) = (3^(-k))/2*sqrt(j/9+n-nFix-j);

       %e_i    = [zeros(i-1,1);1;zeros(n-i,1)];
       e_i    = i;
       Split(i,r) = Split(i,r)+1;

       rightchild = 1; % flag if there will be a right/left child or not,
       leftchild  = 1; % used when splitting along an integer dimension.

       % ******* LEFT NEW POINT ********
       Split_left = Split(:,r);
       if isMIP
          I_L_left = I_L(:,r);
          I_U_left = I_U(:,r);
       end   
       if I_logic(i) % We shall split along an integer dimension
          I_i   = find(IntVars==i);
          aa = I_L(I_i,r);
          bb = I_U(I_i,r);
          delta = floor( (bb-aa+1)/3 );
          c_left = C(:,r);
          if delta >= 1
             I_L_left(I_i) = aa;
             I_U_left(I_i) = aa+delta-1;
             I_L(I_i,r) = aa + delta;
             I_U(I_i,r) = bb - delta;
          elseif delta == 0 % Now there will be only 1 child. Left or right?
             %parent_center = (c_left(i)-x_L(i))./max(x_D(i),1E-20);
             parent_center = x_L(i) + c_left(i)*x_D(i);
             if abs(aa-parent_center) < 1E-4 % if aa==parent_center
                leftchild  = 0;
                rightchild = 0;
                I_L_left(I_i) = bb;
                I_U_left(I_i) = bb;
                I_L(I_i,r) = aa;
                I_U(I_i,r) = aa;
             else
                rightchild=0;
                I_L_left(I_i) = aa;
                I_U_left(I_i) = aa;
                I_L(I_i,r) = bb;
                I_U(I_i,r) = bb;
             end
          else
             rightchild = 0;
             disp('Error in glcSolve, this should not happen');
             return;
          end
          x_i_mid = floor((I_U_left(I_i)+I_L_left(I_i))/2);
          c_left(i) = (x_i_mid-x_L(i))/max(x_D(i),1E-20);
          if I_L_left(I_i)==I_U_left(I_i)
             Split_left(i) = Inf;
if DEBUG
'left'
Iter
e_i
i
r
disp(' ')
pause
end
          end
          if I_L(I_i,r)==I_U(I_i,r)
             Split(i,r) = Inf;
if DEBUG
'node'
Iter
e_i
I_i
r
disp(' ')
pause
end
          end
         %         if length(I)==n
          %if all(~isfinite(Split(:,r)))
          %   ignoreidx = [ignoreidx r];
          %end
          %if all(~isfinite(Split_left))
          %   ignoreidx = [ignoreidx length(F)+1];
          %end
         %         end
       else
          T(r) = T(r)+1;
          delta  = 3^(-Split(i,r));
          if delta < xTol
             Split(i,r)      = inf;
             Split_left(i) = inf;
          end
          %c_left = C(:,r) - delta*e_i; % Centerpoint for new left rectangle
          % Centerpoint for new left rectangle
          c_left = C(:,r);
          c_left(e_i) = c_left(e_i) - delta; 
       end
       if all(isinf(Split(:,r)))
          ignoreidx = [ignoreidx r];
       end
       if all(isinf(Split_left))
          ignoreidx = [ignoreidx length(F)+1];
       end
      
       %if DEBUG & C(i,r)==c_left(i)
       %   fprintf('\n\n Error in glcSolve, C(i,r)==c_left(i) !!!!!! \n');
       %   i   
       %   delta   
       %   Iter       
       %   C(:,r)
       %   c_left
       %   pause
       %end


       % Transform c_left to original search space
       % x_left = x_L + c_left.*x_D;  
       x_left = tomsol(9,x_L, c_left,x_D); 
       f_left = nlp_f(x_left, Prob, varargin{:});   % Function value at x_left
       nFunc  = nFunc+1;
       if LinI + LinE > 0 
          Ax_left = A*x_left;
       else
          Ax_left = [];
       end   
       if NLinI + NLinE > 0
          cx_left = nlp_c(x_left, Prob, varargin{:});
          cx_left = cx_left(:);
       else
          cx_left = [];
       end
       if m > 0
          gx_left    = [-Ax_left(b_L_idx);Ax_left(b_U_idx); ...
                        -cx_left(c_L_idx);cx_left(c_U_idx); ... 
                         Ax_left(leq);cx_left(nleq)]; 
          tmpL = gx_left - g_U;
          if mE > 0
             rE = sum(beta.*abs(tmpL(mI+1:m)));
          else
             rE = 0;
          end

          G(:,nFunc+nFuncOld) = tmpL; % Subtract g_U to get g(x)<=0

          %if all( gx_left(1:mI) < g_U(1:mI) + cTol ) % New point feasible ?

          if all( tmpL(1:mI) <= g_T(1:mI) ) % New point feasible ?
             if feasible < 0
                glcfMin  = f_left + rE; % first feasible point
                feasible = 1;
                fMinIdx  = nFunc;
                xBest    = x_left;
                fMinEQ    = rE;
             elseif f_left + rE < glcfMin % Update glcfMin
                if glcfMin - (f_left + rE) < fEqual
                   % Close point, add to set of minimum points
                   fMinIdx = [nFunc,fMinIdx];
                   xBest   = [x_left,xBest];
                   fMinEQ   = [rE,fMinEQ];
                else
                   glcfMin = f_left + rE;
                   fMinIdx = nFunc;
                   xBest   = x_left;
                   fMinEQ   = [rE,fMinEQ];
                end
             elseif f_left + rE < glcfMin + fEqual % Update fMinIdx
                % Close point, add to set of minimum points
                fMinIdx = [fMinIdx,nFunc];
                xBest   = [xBest,x_left];
                fMinEQ   = [fMinEQ,rE];
             end
          elseif all( tmpL(1:LinI) <= g_T(1:LinI) ) 
             % Linear constraints feasible ?
             feasible = 0;
          end
       else
          if f_left < glcfMin % Update glcfMin
             if glcfMin - (f_left + rE) < fEqual
                % Close point, add to set of minimum points
                fMinIdx = [nFunc,fMinIdx];
                xBest   = [x_left,xBest];
                fMinEQ   = [rE,fMinEQ];
             else
                glcfMin = f_left;
                fMinIdx = nFunc;
                xBest   = x_left;
                fMinEQ   = rE;
             end
          elseif f_left < glcfMin + fEqual-6 % Update fMinIdx
             % Close point, add to set of minimum points
             fMinIdx = [fMinIdx,nFunc];
             xBest   = [xBest,x_left];
             fMinEQ   = [fMinEQ,rE];
          end      
       end
       % Update D(r)
       if nCont > 0
          j = mod(T(r),nCont);
          k = (T(r)-j)/(nCont);
       end
       if isMIP
          if nCont > 0
             z = (3^(-2*k))/4*(j/9+nCont-j);
          else
             z = 0;
          end
          % Length for mid point rectangle
          d = ceil(0.5*(I_U(:,r) - I_L(:,r)))./max(1,x_D(IntVars));
          D(r) = sqrt(z + sum(d.^2));
          %d = (I_U(:,r) - I_L(:,r));
          %d(d==1) = 2;    % Correct length 1 integer variables to length 1
          %D(r) = sqrt(z + sum(d.^2/4));

          % Length for left child rectangle
          d = ceil(0.5*(I_U_left - I_L_left))./max(1,x_D(IntVars));
          Dl = sqrt(z + sum(d.^2));
          %d = (I_U_left - I_L_left);
          %d(d==1) = 2;    % Correct length 1 integer variables to length 1
          %Dl = sqrt(z + sum(d.^2/4));
          D    = [D Dl];
          I_L = [I_L I_L_left];
          I_U = [I_U I_U_left];
       else
          D(r) = (3^(-k))/2*sqrt(j/9+nCont-j);
          D     = [D D(r)];
       end

       C     = [C c_left];
       F     = [F f_left];
       Split = [Split Split_left];
       %T     = [T T(r)];
       T(nFunc+nFuncOld) = T(r);
      
       if feasible < 1, fMin=Inf; else fMin=glcfMin; end
      
       %if convflag == 0
       %   convflag = isClose(fGoal,fMin,fTol,nFunc,Iter,EpsGlob,PriLev);
       %end

       if rightchild
          % ******* RIGHT NEW POINT ********
          Split_right = Split(:,r);
          if isMIP
             I_L_right = I_L(:,r);
             I_U_right = I_U(:,r);
          end   
          i = e_i; % Just for safety, reset i
          if I_logic(i) % We shall split along an integer dimension
             c_right = C(:,r);
             I_L_right(I_i) = bb-delta+1;
             I_U_right(I_i) = bb;
            
             x_i_mid = floor((I_U_right(I_i)+I_L_right(I_i))/2);
             %c_right(i) = (x_i_mid-x_L(i))/max(x_D(i),1E-20);
             % Not possible to get a rightchild if x_D(i) == 0
             c_right(i) = (x_i_mid-x_L(i))/x_D(i);
             %HKHif I_L_left(I_i)==I_U_left(I_i)
             if I_L_right(I_i)==I_U_right(I_i)
                Split_right(i) = Inf;
if DEBUG
'right'
Iter
e_i
I_i
r
disp(' ')
pause
end
             end
         else
            %c_right = C(:,r) + delta*e_i; % Centerpoint for new right rectangle
            % Centerpoint for new right rectangle
            c_right      = C(:,r);
            c_right(e_i) = c_right(e_i) + delta;
         end
         if all(isinf(Split_right))
            ignoreidx = [ignoreidx length(F)+1];
         end
         % Transform c_right to original search space
         %x_right = x_L + c_right.*x_D;  
         x_right = tomsol(9,x_L, c_right,x_D); 
         f_right = nlp_f(x_right, Prob, varargin{:}); % f(x)-value at x_right
         nFunc  = nFunc+1;
         if LinI > 0 | ~isempty(leq)
            Ax_right = A*x_right;
         else
            Ax_right = [];
         end
         if NLinI > 0 | ~isempty(nleq)
            cx_right = nlp_c(x_right, Prob, varargin{:});
            cx_right = cx_right(:);
         else
            cx_right = [];
         end
         if m > 0
            gx_right    = [-Ax_right(b_L_idx);Ax_right(b_U_idx); ...
                           -cx_right(c_L_idx);cx_right(c_U_idx); ... 
                            Ax_right(leq);cx_right(nleq)]; 

            tmpR = gx_right - g_U;

            if mE > 0
               rE = sum(beta.*abs(tmpR(mI+1:m)));
            else
               rE = 0;
            end

            G(:,nFunc+nFuncOld) = tmpR; % Subtract g_U to get g(x)<=0

            %if all( gx_right(1:mI) < g_U(1:mI) + cTol ) % New point feasible ?

            if all( tmpR(1:mI) <= g_T(1:mI) ) % New point feasible ?
               if feasible < 1
                  glcfMin  = f_right + rE; % first feasible point
                  feasible = 1;
                  fMinIdx  = nFunc;
                  xBest    = x_right;
                  fMinEQ    = rE;
               elseif f_right + rE < glcfMin % Update glcfMin
                  if glcfMin - (f_right + rE) < fEqual
                     % Close point, add to set of minimum points
                     fMinIdx = [nFunc,fMinIdx];
                     xBest   = [x_right,xBest];
                     fMinEQ   = [rE,fMinEQ];
                  else
                     glcfMin = f_right + rE;
                     fMinIdx = nFunc;
                     xBest   = x_right;
                     fMinEQ    = rE;
                  end
               elseif f_right + rE < glcfMin + fEqual % Update fMinIdx
                  % Close point, add to set of minimum points
                  fMinIdx = [fMinIdx,nFunc];
                  xBest   = [xBest,x_right];
                  fMinEQ   = [fMinEQ,rE];
               end
            elseif all( tmpL(1:LinI) <= g_T(1:LinI) ) 
               % Linear constraints feasible ?
               feasible = 0;
            end
         elseif f_right < glcfMin % Update glcfMin and x_min
            if glcfMin - (f_right + rE) < fEqual
               % Close point, add to set of minimum points
               fMinIdx = [nFunc,fMinIdx];
               xBest   = [x_right,xBest];
               fMinEQ   = [rE,fMinEQ];
            else
               glcfMin = f_right;
               fMinIdx = nFunc;
               xBest   = x_right;
               fMinEQ   = rE;
            end
         elseif f_right < glcfMin + fEqual % Update fMinIdx
            % Close point, add to set of minimum points
            fMinIdx = [fMinIdx,nFunc];
            xBest   = [xBest,x_right];
            fMinEQ   = [fMinEQ,rE];
         end   
         if isMIP
            I_L = [I_L I_L_right];
            I_U = [I_U I_U_right];
            % Length for right child rectangle
            d = ceil(0.5*(I_U_right - I_L_right))./max(1,x_D(IntVars));
            Dr = sqrt(z + sum(d.^2));
            %d = (I_U_right - I_L_right);
            %d(d==1) = 2;    % Correct length 1 integer variables to length 1
            %Dr = sqrt(z + sum(d.^2/4));
            D    = [D Dr];
         else
            D     = [D D(r)];
         end
         C     = [C c_right];
         F     = [F f_right];
         Split = [Split Split_right];
         %T     = [T T(r)];
         T(nFunc+nFuncOld) = T(r);
         
         if feasible < 1, fMin=Inf; else fMin=glcfMin; end

         %if convflag == 0
         %   convflag = isClose(fGoal,fMin,fTol,nFunc,Iter,EpsGlob,PriLev);
         %end

      end % if rightchild
   
      if m > 0
         % UPDATE THE s(j):s
         % Transform C(:,r) to original search space
         %x_mid = x_L + C(:,r).*x_D;  
         x_mid = tomsol(9,x_L, C(:,r),x_D); 
         norm_left  = norm(x_left-x_mid);
         if AvFunc < AvIter
            if norm_left ~= 0 
               s_0 = (AvFunc*s_0 + abs(f_left-F(r))/norm_left)/(AvFunc+1);
               if NLinI > 0
                  ll=LinI+1;
                  s(ll:mI) = (AvFunc*s(ll:mI) + abs(tmpL(ll:mI) ...
                       - G(ll:mI,r))/norm_left)/(AvFunc+1);
               end
               if NLinE > 0
                  ll = mI+LinE+1;
                  s(ll:m) = (AvFunc*s(ll:m) + abs( tmpL(ll:m) - ...
                             G(ll:m,r)) / norm_left)/(AvFunc+1);
                  %s(ll:m) = (AvFunc*s(ll:m) + abs( abs(tmpL(ll:m)) - ...
                  %           abs(G(ll:m,r))) / norm_left)/(AvFunc+1);
               end
               AvFunc = AvFunc + 1;
            end
            if rightchild
               norm_right = norm(x_right-x_mid);
               if norm_right ~= 0 
                  s_0 = (AvFunc*s_0 + abs(f_right-F(r))/norm_right)/(AvFunc+1);
                  if NLinI > 0
                     %s(LinI+1:m) =(AvFunc*s(LinI+1:m)+abs(gx_right(LinI+1:m)...
                     %             -G(LinI+1:m,r))/norm_right)/(AvFunc+1);
                     ll=LinI+1;
                     s(ll:mI) =(AvFunc*s(ll:mI)+abs(tmpR(ll:mI) ...
                          - G(ll:mI,r))/norm_right)/(AvFunc+1);
                  end
                  if NLinE > 0
                     ll = mI+LinE+1;
                     s(ll:m) =(AvFunc*s(ll:m)+abs( tmpR(ll:m) - ...
                               G(ll:m,r))/norm_right)/(AvFunc+1);
                     %s(ll:m) =(AvFunc*s(ll:m)+abs( abs(tmpR(ll:m)) - ...
                     %          abs(G(ll:m,r)))/norm_right)/(AvFunc+1);
                  end
                  AvFunc = AvFunc + 1;
               end
            end
         else
            % Use exponential smoothing
            alfa = 0.99;
            if norm_left ~= 0 
               s_0         = alfa*s_0 + (1-alfa)*abs(f_left-F(r))/norm_left;
               if NLinI > 0
                  ll=LinI+1;
                  s(ll:mI) = alfa*s(ll:mI) + (1-alfa) * abs(...
                     tmpL(ll:mI)-G(ll:mI,r))/norm_left;
               end
               if NLinE > 0
                  ll = mI+LinE+1;
                  s(ll:m) = alfa*s(ll:m) + (1-alfa) * abs(...
                     abs(tmpL(ll:m))-abs(G(ll:m,r)))/norm_left;
               end
               AvFunc = AvFunc + 1;
            end
            if rightchild
               norm_right = norm(x_right-x_mid);
               if norm_right ~= 0 
                  s_0        = alfa*s_0 + (1-alfa)*abs(f_right-F(r))/norm_right;
                  if NLinI > 0
                     ll=LinI+1;
                     s(ll:mI)= alfa*s(ll:mI) + (1-alfa) * abs(...
                     tmpR(ll:mI)-G(ll:mI,r))/norm_right;
                  end
                  if NLinE > 0
                     ll = mI+LinE+1;
                     s(ll:m)= alfa*s(ll:m) + (1-alfa) * abs( ...
                        abs(tmpR(ll:m))-abs(G(ll:m,r)))/norm_right;
                  end
                  AvFunc = AvFunc + 1;
               end
            end
         end
         if PriLev > 2 
            fprintf('AvF %d s_0 %f',AvFunc,s_0)
            if LinI > 0
               fprintf(' Lin ')
               fprintf(' %f',s(1:LinI))
            end
            if NLinI > 0
               fprintf(' NonLin ')
               fprintf(' %f',s(LinI+1:mI))
            end
            if LinE > 0
               fprintf(' LinE')
               fprintf(' %f',s(mI+1:mI+LinE))
            end
            if NLinE > 0
               fprintf(' NonLinE')
               fprintf(' %f',s(mI+LinE+1:m))
            end
            fprintf('\n')
         end
      end
      if PriLev > 1
         fprintf('Iter:%4d Ev:%4d fMin:%16.6f',...
               Iter, nFunc, fMin);
         if rightchild
            fprintf(' %16.6f %16.6f',F(nFunc-1:nFunc));
         else
            fprintf(' %16.6f                 ',F(nFunc));
         end
         fprintf(' Var %2d',e_i);
         fprintf(' D:%8.4g',delta);
         fprintf(' Split %4d',Split(e_i,r));
         if rightchild
            fprintf(' %4d %4d',Split(e_i,nFunc-1:nFunc));
         else
            fprintf(' %4d    ',Split(e_i,nFunc));
         end
         if ~isempty(IntVars)
            fprintf(' IV:')
            for i = 1:length(IntVars)
                fprintf(' %d',x_left(IntVars(i))); 
            end
            if rightchild
               fprintf(' IV-r:')
               for i = 1:length(IntVars)
                   fprintf(' %d',x_right(IntVars(i))); 
               end
            end
         end
         if isinf(fMin)
            fprintf(' Feas %d',feasible);
            if m > 0, fprintf(' c:%20.10f',min(max(G(:,1:nFunc)))); end
         end
         fprintf('\n')
      end
      
   end % for each rectangle in S
   
   %
   %  STEP 5, UPDATE S (this step is handled by the for loop)
   %
   
   %
   %  STEP 6, ITERATE (this step is handled by the for loop)
   %
   
   if isempty(S) % If no feasible integer solution exist
      break;
   else
      
   if PriLev > 0 | (IterPrint > 0 & glcfMin < glcfOld)   
      fprintf('Iter:%6d  Ev:%6d fMin:%20.10f',...
               Iter, nFunc, fMin);
      %fprintf(' D:%8.4g',delta);
      %fprintf(' Var %2d',e_i);
      %fprintf(' Split %3d',Split(e_i,nFunc));
      if mE > 0
         fprintf(' Eq:%12.10f',fMinEQ(1));
      end
      if ~isempty(IntVars) 
         fprintf(' IV:')
         for i = 1:length(IntVars)
             fprintf(' %d',xBest(IntVars(i),1)); 
         end
      end
      fprintf(' ');
      fprintf(' %f',xBest(find(~I_logic),1)); 
      if isinf(fMin)
         fprintf(' Feas %d',feasible);
         if m > 0, fprintf(' c:%20.10f',min(max(G(:,1:nFunc)))); end
      end
      fprintf('\n');
      glcfOld = glcfMin;
   end

   if convflag == 0
      convflag = isClose(fGoal,fMin,fTol,nFunc,Iter,EpsGlob,PriLev);
   end
   
   end
   
end % while 1  (main loop)
   
fprintf('\n\n');


% SAVE RESULTS
Result = ResultDef(Prob);
Result.Solver = 'glcSolve';
Result.SolverAlgorithm = 'Constrained DIRECT (Jones 2000)';

% Maximal size of possible triangles
tmpset1=ones(1,length(F));
tmpset1(ignoreidx) = 0; % not consider rects being fathomed
Result.maxTri = 3.^-(min(min(Split(:,find(tmpset1))))); 

if feasible > 0
   Result.f_k  = glcfMin;    % Best function value
else
   Result.f_k  = Inf;      % No feasible point found
end
Result.Iter = Iter;     % Number of iterations

%CC = [];
%for i = 1:length(F) % Transform to original coordinates
%   CC = [CC x_L+C(:,i).*x_D];
%end

% For restart

% Save C in transformed form

if m > 0, G = G(:,1:nFuncOld+nFunc); end

save('glcSave.mat','Name','C','F','T','D','G','I_L','I_U','s_0','s','t',...
     'ignoreidx','feasible','Split','glcfMin','fMinIdx','fMinEQ');


% Find all points i with F(i)=glcfMin
if feasible > 0
   Result.minPnts = fMinIdx;
   Result.x_k     = xBest;

   %idx = find(F==glcfMin);
   %if m > 0 % If there are constraints, pick out the feasible points in idx.
   %   if mE == 0
   %      idx2 = all(G(:,idx) < cTol ); % if feasible
   %   elseif mI == 0
   %      idx2 = all(abs(G(:,idx)) < cTol ); % if feasible
   %   else
   %      idx2 = all(G(1:mI,idx) < cTol) & all(abs(G(mI+1:m,idx)) < cTol);
   %   end
   %   idx = idx(idx2);
   %end
   %% All points i with F(i)=glcfMin, and feasible
   %Result.x_k = tomsol(9,x_L,C(:,idx),x_D);    
else
   Result.minPnts=[];
end

if NLinI + NLinE > 0
   c_k=[];
   for i=1:length(Result.minPnts)
      cxx=nlp_c(Result.x_k(:,i), Prob, varargin{:});
      c_k = [c_k  cxx(:)];
   end
   Result.c_k = c_k; % Constraint value at x_k;
end

m=nFunc;
if feasible > 0
   Result.ExitFlag = 0;
   Result.ExitText = ['Tried ' num2str(m+nFuncOld) ' function values in total'];
else
   Result.ExitFlag = 7;
   Result.ExitText = ['Tried ' num2str(m+nFuncOld) ...
                      ' f(x) values in total, no feasible point'];
end
Result.FuncEv=nFunc;

Result=endSolve(Prob,Result);

%if 0  & 1%n==2
%   hold off
%   if 0
%      u = linspace(Prob.x_L(1),Prob.x_U(1),50);
%      v = linspace(Prob.x_L(2),Prob.x_U(2),50);
%      for j=1:length(u)
%        for i=1:length(v)
%           x=[u(j);v(i)];
%           ff(i,j)=nlp_c(x,Prob, varargin{:});
%        end
%      end
%%      contour(u,v,ff,[0.75 1 2 3 4 5 6 7]);
%      hold on
%   end
%%   hold off
%   plot(Result.GLOBAL.C(1,:),Result.GLOBAL.C(2,:),'*');
%   if 0
%      pause 
%      plot(CC(1,1:floor(size(Result.GLOBAL.C,2)/2)),Result.GLOBAL.C(2,1:floor(size(CC,2)/2)),'*');
%      hold on
%      plot(CC(1,floor(size(CC,2)/2)+1:size(CC,2)),CC(2,floor(size(CC,2)/2)+1:size(CC,2)),'o');
%   end
%   grid on
%   zoom on
%end  

function convflag = isClose(fGoal,f,fTol,nFunc,Iter,EpsGlob,PriLev)

convflag = 0;
if isempty(fGoal), return, end
if isinf(fGoal),   return, end

if f <= fGoal
   convflag = 1;
elseif fGoal == 0
   if abs(f-fGoal) < fTol
      convflag = 2;
   end
elseif abs(f-fGoal) <= abs(fGoal) * fTol
   convflag = 3;
end

if convflag > 0 & PriLev >= 0 
   if convflag == 1
      fprintf('\n\nFunction value %f is less than fGoal %f \n',f,fGoal);
   elseif convflag == 2
      fprintf('\n\nError in function value %f is ',f);
      fprintf('%f <= fTol %f\n',abs(f-fGoal),fTol);
   elseif convflag == 3
      fprintf('\n\nRelative error in function value %f is ',f);
      fprintf('%f <= fTol %f\n',abs(f-fGoal)/abs(fGoal),fTol);
   end
   fprintf('Number of function evaluations:  %d\n',nFunc);
   fprintf('Number of iterations:            %d\n',Iter);
   fprintf('Epsilon:                         %f\n',EpsGlob);
end

% MODIFICATION LOG:
%
% 990408  mbk  Modified to make restart possible.
% 990416  mbk  Output if PriLev > 0.
% 990416  mbk  Small changes in comments.
% 000830  hkh  Some speedups
% 000928  hkh  Revision for v3.0
% 001011  hkh  Include Name in glcSave.mat for safety
% 010330  hkh  Check for if ~isempty(tmpidx), row 449, to avoid 6.0 warning
% 010416  hkh  Define isClose, and target value convergence, tomsol speedups 
% 010416  hkh  Bug in speedup of the idxspeed computation
% 011031  hkh  Improve comments. Add maxTri calculation.
% 011110  hkh  Fixed errors in comments
% 011212  hkh  I(tmpidx) should be IL(tmpidx), in error printout
% 020110  hkh  Change name f_min to glcfMin, f_minn to fMin
% 020309  hkh  Clean up code. New way of computing s_0, rateOfChange
% 020311  hkh  Fixed old bug computing rateOfChange. Weight linear/nonlinear
% 020313  hkh  Major revision, speedup of bottle necks.
% 020313  hkh  Change logic to check both slope and constant function
% 020325  hkh  Major revision
% 020325  hkh  Add LinWeight (linear /nonlinear balance) parameter
% 020325  hkh  Add fEqual parameter, tolerance when points are equal
% 020420  hkh  Incorrect check, checking left. not right int triangle
% 040415  hkh  Change iniSolve call
% 060814  med  FUNCS used for callbacks instead
