% glbDirect implements a modified version of the algorithm DIRECT.
%
% Reference: D. R. Jones, C. D. Perttunen and B. E. Stuckman,
% Lipschitzian Optimization Without the Lipschitz Constant,
% JOTA  Vol. 79, No. 1, October 1993.
%
% glbDirect solves box-bounced global optimization problems of the form:
%
%    min   f(x)
%     x
%    s/t   x_L <= x <= x_U
%
% Calling syntax:
%
% function Result = glbDirectTL(Prob)
%
% See below "USAGE" on how to create the Prob structure and do the call
% Read the part: "IMPORTANT NOTE ABOUT THE DIRECT ALGORITHM"
%
% INPUT PARAMETERS
%
% Prob        Structure, where the following variables are used:
%   Name      Name of the problem. Used for security if doing warm start.
%   FUNCS.f   The routine to compute the function, given as a string.
%             A call to tomFiles.m or glcAssign.m is setting this field.
%   x_L       Lower bounds for each element in x.
%             Any lower bounds that are inf are changed to -10000
%   x_U       Upper bounds for each element in x.
%             Any upper bounds that are inf are changed to  10000
%   PriLevOpt Print Level
%                0 = Silent
%                1 = Errors
%                2 = Termination message and warm start info
%                3 = Option summary
%   WarmStart If true, >0, glbDirect reads the output from the last
%             run from Prob.glbDirect.WarmStartInto if it
%             exists. If it doesn't exist, glbDirect attempts to
%             open and read warm start data from mat-file
%             glbDirectSave.mat. glbDirect uses this warm start
%             information to continue from the last run.
% optParam    Structure in Prob, Prob.optParam
%             Defines optimization parameters. Fields used:
%  IterPrint  Print iteration log every IterPrint iteration. Set to
%             0 for no iteration log. PriLev must be set to at
%             least 1 to have iteration log to be printed.
%  MaxIter    Maximal number of iterations, default 200.
%  MaxFunc    Maximal number of function evaluations, default 10000
%  EpsGlob    Global/local weight parameter, default 1E-4.
%  fGoal      Goal for function value, if empty not used
%  eps_f      Relative accuracy fTol for function value, fTol == eps_f
%             Stop if abs(f-fGoal) <= abs(fGoal) * fTol , if fGoal~=0
%             Stop if abs(f-fGoal) <= fTol , if fGoal==0
%  eps_x      Convergence tolerance in x. All possible rectangles are
%             less than this tolerance. Scaled automatically.
%             See the output field maxTri.
%
% glbDirect   Structure in Prob, Prob.glbDirect. Solver specific.
%  options    Structure with options. These options have precedence
%             over all other options in the Prob struct. Available
%             options are:
%
%               PRILEV    Equivalent to Prob.PrilevOpt.
%                         Default: 0
%
%               MAXFUNC   Eq. to Prob.optParam.MaxFunc.
%                         Default: 10000
%
%               MAXITER   Eq. to Prob.optParam.MaxIter.
%                         Default: 200
%
%               PARALLEL  Set to 1 in order to have glbDirect to
%                         call Prob.FUNCS.f with a matrix x of points
%                         to let the user function compute function
%                         values in parallel.
%                         Default: 0
%
%               WARMSTART Eq. to Prob.WarmStart.
%                         Default: 0
%
%               ITERPRINT Eq. to Prob.optParam.IterPrint.
%                         Default: 0
%
%               FUNTOL    Eq. to Prob.optParam.eps_f.
%                         Default: 1e-2
%
%               VARTOL    Eq. to Prob.optParam.eps_x.
%                         Default: 1e-13
%
%               GLWEIGHT  Eq. to Prob.optParam.EpsGlob.
%                         Default: 1e-4
%
%  WarmStartInfo
%             Structure with WarmStartInfo. Use WarmDefGLOBAL.m to
%             define it.
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
%               or Relative error in function value f is less than
%               fTol, i.e. abs(f-fGoal)/abs(fGoal) <= fTol
%           3 = Maximum number of iterations done
%           4 = Maximum number of function evaluations done
%           91= Numerical trouble, did not find element in list
%           92= Numerical trouble, No rectangle to work on
%           99= Other error, see ExitFlag
%
%  glbDirect  Substructure for glbDirect specific result data.
%    nextIterFunc
%             If optimization algorithm was stopped because of
%             maximum number of function evaluations reached, this is
%             the number of function evaluations required to complete
%             the next iteration.
%    maxTri   Maximum size of any triangles.
%    WarmStartInfo
%             Structure containing warm start data. Could be used
%             to continue optimization where glbDirect stopped.
%
% To make a warm start possible, glbDirect saves the following information in
% the structure Result.glbDirect.WarmStartInfo and file
% glbDirectSave.mat (for internal solver use only):
%   points    Matrix with all rectangle centerpoints, in [0,1]-space.
%   dRect     Vector with distances from centerpoint to the vertices.
%   fPoints   Vector with function values.
%   nIter     Number of iterations.
%   lRect     Matrix with all rectangle side lengths in each dimension.
%   Name      Name of the problem. Used for security if doing warm start.
%   dMin      Row vector of minimum function value for each distance.
%   ds        Row vector of all different distances, sorted.
%   glbfMin   Best function value found at a feasible point.
%   iMin      The index in D which has lowest function value, i.e. the
%             rectangle which minimizes (F - glbfMin + E)./D where
%             E = max(EpsGlob*abs(glbfMin),1E-8)
%   ign       Rectangles to be ignored in the rect. selection procedure.
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
%      % Another way of setting options is the following. However,
%      % these will ONLY apply to glbDirect. These are solver
%      % specific that is:
%      Prob.glbDirect.options.FUNTOL=1e-5;
%
% -- Calling glbDirect
%
% Direct solver call:
%      Result = glbDirectTL(Prob);
%      PrintResult(Result);
%
% Driver call, including printing with level 2:
%      Result = tomRun('glbDirect',Prob,2);
%
% -- Defining an objective function
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
% glbDirect has an option to pass a matrix x of several points for
% making function evaluations in parallel. For heavy calculations,
% a user could divide the points to several threads/computers in
% order to speed up the function evaluation. See option: PARALLEL.
%
% -- Warm start
%
% To make a restart, just set the restart flag, and call glbDirect
% once again to call reading glbDirectSave.mat file:
%
%      Prob.WarmStart = 1;
%      Result = tomRun('glbDirect', Prob, 2);
%
% Another warm start (with same MaxFunc) is made by just calling tomRun again
%      Result = tomRun('glbDirect', Prob, 2);
%
% To make a restart from the warm start information in the Result
% structure, make a call to WarmDefGLOBAL before calling glbDirect:
%
%      Prob = WarmDefGLOBAL('glbDirect', Prob, Result);
%      Result = tomRun('glbDirect', Prob, 2);
% where Result is the result structure returned by the previous run.
% To make another warm start (with same MaxFunc), repeat the two lines:
%      Prob = WarmDefGLOBAL('glbDirect', Prob, Result);
%      Result = tomRun('glbDirect', Prob, 2);
%
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
%
% IMPORTANT NOTE ABOUT THIS IMPLEMENTATION OF THE DIRECT ALGORITHM:
%
% glbDirect never stops optimzation in the middle of an
% iteration. This is to make sure warm start information will be
% complete when returning. This often results in the behavior
% where the solver returns with inform = 2 and the number of
% function evaluations (nFunc) has not reached maxFunc. The reason
% is that the next iteration would require more function evaluations
% than the remaining number of function evaluations. When glbDirect
% returns with inform = 2, the number of function evaluations
% requried for the next iteration is stored in the output variable:
% Result.glbDirect.nextIterFunc.
%
% A special case is when glbDirect returns with nIter unchanged
% (equal to 0 for cold start, and equal to the
% Result.glbDirect.WarmStartInfo.nIter of the previous run. Then no
% new iterations were done. Set Prob.optParam.MaxFunc to at least
% Result.glbDirect.nextIterFunc.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Mar 9, 2005.    Last modified Oct 2, 2009.

function Result = glbDirectTL(Prob)

%#function lp_f lp_g lp_H

if nargin < 1
   error('glbDirect needs input structure Prob');
end

solvType=checkType('glb');

Prob=ProbCheck(Prob,'glbDirect',solvType);

Prob = iniSolve(Prob,solvType,0,0);

%if isempty(Prob.x_L) | isempty(Prob.x_U)
%   disp('glbDirect requires both lower and upper variable bounds');
%   Result.ExitFlag = 1;
%   Result.ExitText = 'glbDirect requires both lower and upper variable bounds';
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

% Safeguard
if MaxIter < 0
   MaxIter = [];
end
if MaxFunc < 0
   MaxFunc = [];
end

Result                 = ResultDef(Prob);
Result.Solver          = 'glbDirect';
Result.SolverAlgorithm = 'DIRECT - Fortran Implementation';

% This code is never executed, as the bounds are changed to -10000
% and 10000 for all infinite elements.
%if any(isinf(x_L)) | any(isinf(x_U))
%   disp('glbDirect solves box-bounded problems.');
%   disp('Found some bound to be Inf');
%   Result.ExitFlag = 2;
%   Result.Inform   = 99;
%   Result.ExitText = str2mat('glbDirect solves box-bounded problems' ...
%                            ,'Found some bound to be Inf');
%   Result=endSolve(Prob,Result);
%   return
%end

if ~isempty(fGoal)
  if isinf(fGoal)
    fGoal = [];
  end
end

fcn = Prob.FUNCS.f;
if xnargin(fcn) > 1
   probPtr = Prob;
else
   probPtr = [];
end

GLBDIRECTSTRUCT = DefPar(Prob, 'glbDirect', []);
options = DefPar(GLBDIRECTSTRUCT, 'options', []);

if Prob.WarmStart == 1
   % Restart with values from previous run.

   Name1 = deblank(Prob.Name);  % Problem name
   
   WarmStartInfo = DefPar(GLBDIRECTSTRUCT, 'WarmStartInfo', []);
   
   % If the WarmStart structure exists, then use that warm start data.
   if ~isempty(WarmStartInfo)
     Name = WarmStartInfo.Name;
     if strcmp(Name1,Name)
       if(PriLev >= 2)
         fprintf(['Restart with %d sampled points from last ' ...
                  'runs.\nUsing warm start data from structure: ' ...
                  'Prob.glbDirect.WarmStartInfo.\n'],...
                  WarmStartInfo.nFunc);
       end
     else
       Prob.WarmStart = 0;
       if PriLev >= -1000
         fprintf('Previous run was with Problem %s\n',Name);
         fprintf('This run is with Problem %s\n',Name1);
         fprintf(['Impossible to do restart using data from ' ...
                  'structure: Prob.glbDirect.WarmStartInfo.\n']);
       end
     end

   % If no WarmStart struture exists, look for a file.
   else
     filename = 'glbDirectSave.mat';
     
     if(exist(filename,'file')~=2)
       fprintf(['Couldn''t find warm start info structure nor warm ' ...
                'start file in path.\nNo warm start possible.\n']);
       Prob.WarmStart = 0;
     else
       load(filename,'Name');

       if strcmp(Name1,Name)
         WarmStartInfo = load(filename,'nFunc', 'nIter', 'glbfMin', ...
                              'points', 'fPoints', 'dRect', 'lRect', ...
                              'ign', 'ds', 'dMin', 'iMin');       
         if PriLev >= 2
           fprintf(['Restart with %d sampled points from last ' ...
                    'runs.\nUsing warm start data from file: %s.\n'], ...
                   WarmStartInfo.nFunc, filename);
         end
         
       else
         Prob.WarmStart = 0;
         if PriLev >= -1000
           fprintf('Previous run was with Problem %s\n',Name);
           fprintf('This run is with Problem %s\n',Name1);
           fprintf(['Impossible to do restart using data from file: ' ...
                    '%s.\n'], filename);
         end
       end
     end
   end
end

if Prob.WarmStart == 1 
  Iter0  = WarmStartInfo.nIter;
  nFunc0 = WarmStartInfo.nFunc;
else
  Iter0  = 0;
  nFunc0 = 0;
  WarmStartInfo = [];
end

% Set options.

if ~isempty(PriLev)
  options.PRILEV = DefPar(options, 'PRILEV', PriLev);
end

if ~isempty(MaxFunc)
  options.MAXFUNC = DefPar(options, 'MAXFUNC', MaxFunc);
end

if ~isempty(MaxIter)
  options.MAXITER = DefPar(options, 'MAXITER', MaxIter);
end

if ~isempty(Prob.WarmStart)
  options.WARMSTART = DefPar(options, 'WARMSTART', double(Prob.WarmStart > 0));
end

if ~isempty(IterPrint)
  options.ITERPRINT = DefPar(options, 'ITERPRINT', IterPrint);
end

if ~isempty(fTol)
  options.FUNTOL = DefPar(options, 'FUNTOL', fTol);
end

if ~isempty(xTol)
  options.VARTOL = DefPar(options, 'VARTOL', xTol);
end

if ~isempty(EpsGlob)
  options.GLWEIGHT = DefPar(options, 'GLWEIGHT', EpsGlob);
end

Name = deblank(Prob.Name);  % Problem name

[convFlag, xMin, glbfMin, nFunc, Iter, maxTri, nextIterFunc, WarmStartInfoOut] ...
    = glbDirect(fcn, x_L, x_U, fGoal, options, WarmStartInfo, probPtr);

% Add name to WarmStartInfo structure. It is TL-specific.
WarmStartInfoOut.Name = Name;
Result.glbDirect.WarmStartInfo = WarmStartInfoOut;

SaveWarmStartFile(WarmStartInfoOut);

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
 otherwise,
  ExitFlag  = 4;
  Inform    = convFlag;
  Result.ExitText = 'Internal error. Contact tomlab@tomopt.com.';
end

Result.glbDirect.maxTri = maxTri;       % Maximal size of possible
                                        % triangles
Result.maxTri           = maxTri;
Result.glbDirect.nextIterFunc = nextIterFunc;
Result.x_k              = xMin;         % Best point(s) found
Result.f_k              = glbfMin;      % Best function value
Result.Iter             = Iter-Iter0;  % Number of iterations this run
Result.FuncEv           = nFunc-nFunc0; % Number of function evaluations this run
Result.ExitFlag         = ExitFlag;
Result.Inform           = Inform;

% Must set global n_f to see number of function evaluations in the GUI

global n_f
n_f = nFunc-nFunc0;

Result          = endSolve(Prob,Result);

function SaveWarmStartFile(WarmStartInfo)

%     .nFunc
%     .nIter
%     .glbfMin
%     .points
%     .fPoints
%     .dRect
%     .lRect
%     .ign
%     .ds
%     .dMin
%     .iMin
%     plus .Name

filename = 'glbDirectSave.mat';

nFunc   = WarmStartInfo.nFunc;
nIter   = WarmStartInfo.nIter;
glbfMin = WarmStartInfo.glbfMin;
points  = WarmStartInfo.points;
fPoints = WarmStartInfo.fPoints;
dRect   = WarmStartInfo.dRect;
lRect   = WarmStartInfo.lRect;
ign     = WarmStartInfo.ign;
ds      = WarmStartInfo.ds;
dMin    = WarmStartInfo.dMin;
iMin    = WarmStartInfo.iMin;
Name    = WarmStartInfo.Name;

% Safed with try-catch - a failure here will lose the entire run otherwise
try
   save(filename, 'nFunc', 'nIter', 'glbfMin', 'points', 'fPoints', ...
      'dRect', 'lRect', 'ign', 'ds', 'dMin', 'iMin', 'Name');
catch
   warning('Failed to save warmstart information to glbDirectSave.mat');
   disp(lasterr)
   disp('Warm start information is available in Result.glbDirect.WarmStartInfo');
end

% MODIFICATION LOG
%
% 050309 frhe File written, based on glbFast.m.
% 050310 frhe Help modified.
% 060814 med  FUNCS used for callbacks instead
% 061003 ango try-catch safe of save statement
% 091002 hkh  Use WarmDefGLOBAL to generate warm start info in Prob structure
% 091002 hkh  Remove comments about glbFast
