% glbSolve implements a modified version of the algorithm DIRECT.
% This version, glbFast, is a Fortran MEX version of glbSolve
%
% Reference: D. R. Jones, C. D. Perttunen and B. E. Stuckman,
% Lipschitzian Optimization Without the Lipschitz Constant,
% JOTA  Vol. 79, No. 1, October 1993.
%
% glbFast solves box-bounced global optimization problems of the form:
%
%    min   f(x)
%     x
%    s/t   x_L <= x <= x_U   
%
% Calling syntax:
%
% function Result = glbFast(Prob)
%
% See below "USAGE" on how to create the Prob structure and do the call
% Read the part: "IMPORTANT NOTE ABOUT THE DIRECT ALGORITHM"
%
% INPUT PARAMETERS
%
% Prob        Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start
%   FUNCS.f   The routine to compute the function, given as a string.
%             A call to tomFiles.m or glcAssign.m is setting this field. 
%   x_L       Lower bounds for each element in x.
%             Any lower bounds that are inf are changed to -10000
%   x_U       Upper bounds for each element in x.
%             Any upper bounds that are inf are changed to  10000
%   PriLevOpt Print Level 
%             0 = silent. 1 = Warm Start info 2 = A header printed
%   WarmStart If true, >0, glbFast reads the output from the last run from the
%             mat-file glbFastSave.mat, and continues from the last run.
% optParam    Structure in Prob, Prob.optParam 
%             Defines optimization parameters. Fields used: 
%  IterPrint  Print iteration #, # of evaluated points and best f(x) each iter
%  MaxIter    Maximal number of iterations, default max(5000,n*1000);
%  MaxFunc    Maximal number of function evaluations, default max(10000,n*2000)
%  EpsGlob    Global/local weight parameter, default 1E-4.
%  fGoal      Goal for function value, if empty not used
%  eps_f      Relative accuracy fTol for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol , if fGoal~=0
%             Stop if abs(f-fGoal) <= fTol , if fGoal==0
%  eps_x      Convergence tolerance in x. All possible rectangles are 
%             less than this tolerance (scaled to (0,1) )
%             See the output field maxTri.
%
% OUTPUT PARAMETERS
%
% Result    Structure with results from optimization
%  x_k      Matrix with optimal points as columns.
%  f_k      The best function value found so far
%  Iter     Number of iterations
%  FuncEv   Number of function evaluations
%  maxTri   Maximum size of any triangle
%  ExitText Text string giving ExitFlag and Inform information
%  ExitFlag 0 = Normal termination, max number of iterations /func.evals reached
%           1 = Some bound, lower or upper is missing
%           2 = Some bound is inf, must be finite
%           4 = Numerical trouble determining optimal rectangle, empty set
%               and cannot continue
%  Inform   1 = Function value f is less than fGoal
%           2 = Absolute function value f is less than fTol, only if fGoal = 0
%            or Relative error in function value f is less than fTol, i.e.
%               abs(f-fGoal)/abs(fGoal) <= fTol
%           3 = Maximum number of iterations done
%           4 = Maximum number of function evaluations done
%           91= Numerical trouble, did not find element in list
%           92= Numerical trouble, No rectangle to work on
%           99= Other error, see ExitFlag
%
% To make a warm start possible, glbFast saves the following information in
% the file glbFastSave.mat (for internal solver use only):
%   C         Matrix with all rectangle centerpoints, in [0,1]-space.
%   D         Vector with distances from centerpoint to the vertices.
%   E         Computed tolerance in rectangle selection
%   F         Vector with function values.
%   Iter      Number of iterations
%   L         Matrix with all rectangle side lengths in each dimension.
%   Name      Name of the problem. Used for security if doing warm start
%   dMin      Row vector of minimum function value for each distance 
%   ds        Row vector of all different distances, sorted.
%   glbfMin   Best function value found at a feasible point.
%   iMin      The index in D which has lowest function value, i.e. the
%             rectangle which minimizes (F - glbfMin + E)./D where
%             E = max(EpsGlob*abs(glbfMin),1E-8)
%   ignoreIdx Rectangles to be ignored in the rect. selection proceedure.
%
% USAGE:
%
% Let the name of the problem be "GLBF Test"
% The function GLBF is best written as
%     function f = GLBF(x, Prob) 
% Then any information, say b and C is easily sent to GLBF in the
% Prob structure by the call (assume bounds x_L and x_U are initialized) 
%
%      Prob   = glcAssign('GLBF',x_L,x_U,'GLBF Test');
%      Prob.user.b = b; Prob.user.C=C;        % example of extra user data
%
%      % Default values are now set for PriLevOpt, and structure optParam
%      % To change a value, examples are shown on the two next lines
%      Prob.optParam.MaxFunc = 500; % Change max number of function evaluations 
%      Prob.optParam.MaxIter = 100; % Change the number of iterations to 100
%
% Direct solver call:
%      Result = glbFast(Prob);
%      PrintResult(Result);
%
% Driver call, including printing with level 2:
%      Result = tomRun('glbFast',Prob,2);
%           
% The user function GLBF is written
%
%      function f = GLBF(x, Prob) 
%      b = Prob.user.b; C = Prob.user.C;
%      f = "some function of x, b and C"
%
% It is also possible to use the function format
%      function f = GLBF(x) 
%
% but then any additional parameters must be sent as global variables.
%
% To make a restart, just set the restart flag, and call glbFast once again:
%
%      Prob.WarmStart = 1;
%      Result = tomRun('glbFast',Prob,2);
%
% IMPORTANT NOTE ABOUT THE DIRECT ALGORITHM:
%
% The DIRECT algorithm only reaches the variable bounds in the limit.
% Therefore convergence for global optimum where components are on the bounds
% is slow. 
% One remedy is to reduce lower bounds with a tolerance, say 1E-4, and add
% a similar tolerance 1E-4 to the upper bounds that might be reached.
% Another possibility is to fix a variable on its bound by setting the lower
% and upper bounds equal.
%
% Always try to reduce the dimension as much as possible when using the
% DIRECT algorithm, and try to shrink the box, defined by the lower and
% upper bounds, as much as possible.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2006 by Tomlab Optimization Inc., $Release: 5.6.0$
% Written Sep 1, 1998.    Last modified Oct 3, 2006.

function Result = glbFast(Prob)

%#function lp_f lp_g lp_H

if nargin < 1
   error('glbFast needs input structure Prob');
end

solvType=checkType('glb');

Prob=ProbCheck(Prob,'glbFast',solvType);

Prob = iniSolve(Prob,solvType,0,0);

%if isempty(Prob.x_L) | isempty(Prob.x_U)
%   disp('glbFast requires both lower and upper variable bounds');
%   Result.ExitFlag = 1;
%   Result.ExitText = 'glbFast requires both lower and upper variable bounds';
%   Result=endSolve(Prob,Result);
%   return;
%end

% Pick up input parameters from the Prob structure:
x_L = Prob.x_L(:);   % Lower bounds
x_U = Prob.x_U(:);   % Upper bounds

% Check for Inf and set to lower values.
x_L(isinf(x_L)) = -10000;
x_U(isinf(x_U)) =  10000;

PriLev    = Prob.PriLevOpt;          % Print level
MaxIter   = Prob.optParam.MaxIter;   % Number of iterations
MaxFunc   = Prob.optParam.MaxFunc;   % Number of function evaluations
EpsGlob   = Prob.optParam.EpsGlob;   % global/local weight parameter. 
IterPrint = Prob.optParam.IterPrint; % Print short information each iteration
fGoal     = Prob.optParam.fGoal;     % Goal for f(x).
fTol      = Prob.optParam.eps_f;     % Relative tolerance for fGoal
xTol      = Prob.optParam.eps_x;     % Convergence tolerance in x. All possible
                                     % rectangles are less than this tolerance 
                                     % (scaled to (0,1) )
n         = Prob.N;

% Safeguard
if isempty(MaxIter) | MaxIter < 0
   MaxIter = max(5000,n*1000);
end
if isempty(MaxFunc) | MaxFunc < 0
   MaxFunc = max(10000,n*2000);
end

Result                 = ResultDef(Prob);
Result.Solver          = 'glbFast';
Result.SolverAlgorithm = 'DIRECT - Fortran Mex implementation';

if any(isinf(x_L)) | any(isinf(x_U))
   disp('glbFast solves box-bounded problems.');
   disp('Found some bound to be Inf');
   Result.ExitFlag = 2;
   Result.Inform   = 99;
   Result.ExitText =str2mat('glbFast solves box-bounded problems' ...
                           ,'Found some bound to be Inf');
   Result=endSolve(Prob,Result);
   return
end


if isinf(fGoal),   fGoal = -1E300; end
if isempty(fGoal), fGoal = -1E300; end

fcn = Prob.FUNCS.f;
if xnargin(fcn) > 1
   probPtr = Prob;
else
   probPtr = -999;
end

if Prob.WarmStart
   % Restart with values from previous run.

   Name1 = deblank(Prob.Name);  % Problem name
   load('glbFastSave.mat','Name')
   if strcmp(Name1,Name)
      load('glbFastSave.mat','F','Iter')
      nFunc0 = length(F); 
      Iter0 = Iter; 
      clear F
      if PriLev > 0
         fprintf('\n Restart with %d sampled points from last runs\n',nFunc0);
      end
   else
      Prob.WarmStart = 0;
      if PriLev >= -1000
         fprintf('Previous run was with Problem %s\n',Name);
         fprintf('This run is with Problem %s\n',Name1);
         fprintf('Impossible to do restart.\n');
         fprintf('Maybe there exists several files glbFastSave.mat?\n');
      end
      nFunc0 = 0; 
      Iter0 = 0; 
      Name = Name1; 
   end
else
   nFunc0 = 0; 
   Iter0 = 0; 
   Name = deblank(Prob.Name);  % Problem name
end

[glbfMin, xMin, Iter, nFunc, convFlag, maxTri] = ...
   tomsol(16, x_L, x_U, fcn, probPtr, Prob.WarmStart, ...
   MaxFunc, MaxIter, fGoal, fTol, IterPrint, EpsGlob, xTol);

try
   save 'glbFastSave.mat' Name -APPEND
catch
   warning('Failed to append variable Name to glbFastSave.mat');
   disp(lasterr);
end

switch convFlag
 case 3
   ExitFlag = 0;
   Inform   = 3;
   if Prob.WarmStart
      Result.ExitText = ['Max iterations. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(Iter) ' iter.'];
   else
      Result.ExitText = 'Maximal number of iterations reached';
   end
 case 2
   ExitFlag = 0;
   Inform   = 4;
   if Prob.WarmStart
      Result.ExitText = ['Max f(x) evals. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(Iter) ' iter.'];
   else
      Result.ExitText = 'Maximum number of f(x) evaluations reached';
   end
 case 1
   ExitFlag = 0;
   Inform   = 2;
   if Prob.WarmStart
      Result.ExitText = ['Converged to fGoal. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(Iter) ' iter.'];
   else
      Result.ExitText = 'Converged to fGoal';
   end
 case 6
   ExitFlag = 0;
   Inform   = 1;
   if Prob.WarmStart
      Result.ExitText = ['Converged to f(x) < fGoal. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(Iter) ' iter.'];
   else
      Result.ExitText = 'Converged to f(x) < fGoal';
   end
 case 4
   % Numerical trouble, did not find element in list
   ExitFlag = 4;
   Inform   = 91;
   if Prob.WarmStart
      Result.ExitText = ['Numerical trouble - EXIT. ', ...
          'Tried total ' num2str(nFunc) ' f(x), ' num2str(Iter) ' iter.'];
   else
      Result.ExitText = 'Numerical trouble - EXIT!';
   end
 case 5
   % No rectangle to work on
   ExitFlag = 4;
   Inform   = 92;
   Result.ExitText = 'No rectangle to work on - EXIT!';
end

Result.maxTri   = maxTri;       % Maximal size of possible triangles
Result.x_k      = xMin;         % Best point(s) found
Result.f_k      = glbfMin;      % Best function value
Result.Iter     = Iter-Iter0;   % Number of iterations this run
Result.FuncEv   = nFunc-nFunc0; % Number of function evaluations this run
Result.ExitFlag = ExitFlag;
Result.Inform   = Inform;

% Must set global n_f to see number of function evaluations in the GUI

global n_f
n_f = nFunc-nFunc0;

Result          = endSolve(Prob,Result);

%C = C(:,1:m);    % All lengths
%F = F(1:m);      % All function values computed
%D = D(1:m);      % All distances
%L = L(:,1:m);    % All lengths

% Not used
%   t         t(i) is the total # splits along dimension i.

% MODIFICATION LOG
%
% 010601  hkh  See the log for glbSolve
% 010621  jkk  It is all moved to Fortran
% 010715  hkh  Handle more output parameters; create glbSave.mat; test fGoal
% 010815  hkh  Changes for warm start
% 010817  hkh  Avoid calling Tomlab gateway nlp_f, call user
%              routine directly. Return max size of triangles as maxTri
% 011030  hkh  Improve comments. Add Name to mat-file. Check of warm start.
%              Bug fixes for GUI.
% 011112  hkh  Add test on fGoal empty, set as very low then, avoid DLL crash
% 020110  hkh  Change name fMin to glbfMin. Conflict in 5.x with fmin function
% 020506  hkh  Improving comments about DIRECT algorithm
% 040111  hkh  Change call to inisolve
% 040327  hkh  Revised warm start comments
% 040330  hkh  Add Inform value, revised ExitFlag values, test f(x) < fGoal
% 041123  hkh  Change call to tomRun in help
% 050117  med  mlint review
% 051006  hkh  Safe handling of MaxFunc, new default
% 060814  med  FUNCS used for callbacks instead
% 061003  ango try-catch on save statement
