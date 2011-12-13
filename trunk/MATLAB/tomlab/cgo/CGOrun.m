% CGOrun runs the predefined global optimization problems in e.g. glb_prob.m
% for the TOMLAB /CGO solvers rbfSolve and EGO (and other global solvers).
%
% They stop when reaching max number of function evaluations
% or when the the function value goal fGoal is reached with a relative
% precision EpsF.
%
% function Result = CGOrun(P, MaxFunc, Solvers, EpsF, Percent, InitFile)
%
%
% P            Problem number (default 1)
%
% MaxFunc      Maximal number of function evaluations in rbfSolve and ego
%
% Solvers      Select which solvers to run, 0/1 vector of length 9:
%              Solvers(1) = 1  rbfSolve with glcFast as globalSolver (default)
%              Solvers(2) = 1  rbfSolve with glcCluster as globalSolver
%              Solvers(3) = 1  ego with glcFast as globalSolver
%              Solvers(4) = 1  ego with glcCluster as globalSolver
%              Solvers(5) = 1  glbFast (must be unconstrained problem)
%              Solvers(6) = 1  glcFast
%              Solvers(7) = 1  glcCluster
%              Solvers(8) = 1  glbSolve
%              Solvers(9) = 1  glcSolve
% EpsF         Relative precision in achieving fGoal. The goal value fGoal
%              should be set in Prob.f_opt or in Prob.MIP.fGoal.
% 
% Percent      Type of strategy to get the initial sampled values: (default -1)
%
%              See help rbfSolve and help ego for the values of Percent
%
% InitFile     Name of Tomlab init file in the IF format (default 'glb_prob')
%
% Example: Solve Problem 4, Hartman 3, which shows good convergence
%
%       CGOrun(4,60)         Run rbfSolve for 60 function evaluations
%       CGOrun(4,60,[0 0 1]) Run ego for 60 function evaluations
%
% The results are printed on diary file cgoP"z".txt where "z" = Problem number

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Jan 11, 2002.   Last modified Feb 20, 2005.

function Result = CGOrun(P, MaxFunc, Solvers, EpsF, Percent, InitFile)

if nargin < 6
   InitFile = [];
   if nargin < 5
      Percent = [];
      if nargin < 4
         EpsF = [];
         if nargin < 3
            Solvers = [];
            if nargin < 2
               MaxFunc = [];
               if nargin < 1
                  P = [];
end, end, end, end, end, end

if isempty(P),        P = 4; end
if isempty(Solvers),  Solvers = [1; zeros(8,1)]; end
if isempty(EpsF),     EpsF = 1E-5; end
%if isempty(Percent),  Percent = -1; end
if isempty(InitFile), InitFile = 'glb_prob'; end

m = length(Solvers);
if m < 9, Solvers = [Solvers(:); zeros(9-m,1)]; end

Prob = probInit(InitFile,P);     % Initialize using Tomlab IF problem file
Prob.optParam.eps_f     = EpsF;  % Relative precision in reaching fGoal
Prob.optParam.IterPrint = 1;     % Print one line each iteration

diary(['cgoP' num2str(P) '.txt'])

fGoal = Prob.f_opt;

if isempty(fGoal)
   fGoal = Prob.MIP.fGoal;
   fprintf('\nGoal    f(x) %f\n\n',fGoal);
else
   fprintf('\nOptimal f(x) %f\n\n',fGoal);
end
Prob.optParam.fGoal = fGoal;
n = length(Prob.x_L);

if isempty(MaxFunc),  MaxFunc = 20*n; end

Prob.optParam.MaxFunc   = MaxFunc;  % Max function evaluations in rbfSolve/ego

fprintf('\n\nCGOrun: Solve Problem %s, # %d, size %d\n\n',Prob.Name,P,n);

% Special CGO parameters for rbfSolve
Prob.CGO.idea      = 2;  % 2/1
Prob.CGO.rbfType   = 2;  % 2/1
Prob.CGO.SCALE     = 1;
Prob.CGO.PLOT      = 0;
Prob.CGO.REPLACE   = 1;  % 1/0
Prob.CGO.LOCAL     = 0;  % 1/0
Prob.CGO.infStep   = 1;  % 1 = add extra search step with target -inf
Prob.CGO.AddMP     = 0;  % Add midpoint to corner strategy or Gutmann
                         % strategy if AddMP = 1
Prob.CGO.fStarRule = 1;  % 1/2/3. 1 squares weights, 2 = no squaring
                         % 3 = Only using
Prob.CGO.DeltaRule = 1;  % 1 = Skip large f(x) computing interval Delta
                         % 0 = Use all points

% Cycle lengths
if Prob.CGO.idea == 1
   Prob.CGO.N         = 5; % default 5 for fStarRule 1,2. default 1 for rule 3
else
   Prob.CGO.N         = 3; % Always 3
end


Prob.CGO.Percent=Percent;

if Solvers(1) | Solvers(2)
   Prob.PriLevOpt = 0;
   Prob.optParam.IterPrint = 1;     % Print one line each iteration
   Prob.optParam.MaxIter   = 1000;  % Max local iterations in each step
   %Prob.optParam.MaxFunc  = 400;   % Max function evaluations in rbfSolve

   if Solvers(2)
      Prob.CGO.globalSolver = 'glcCluster';
      Prob.GO.maxFunc1 = 401;
      Prob.GO.maxFunc2 = 200;
   else
      Prob.GO.MaxFunc         = 500;   % Max func eval in global search
      Prob.GO.MaxIter         = 500;   % Max iterations in global search
      Prob.CGO.globalSolver = 'glbFast';
      Prob.CGO.globalSolver = 'glcFast';
   end

   % Print info
   CGO=Prob.CGO
   EpsF=Prob.optParam.eps_f

%profile on
   Result = tomRun('rbfSolve',Prob,2);
%profile report

   %Prob.optParam.MaxFunc   = 20;  % Max function evaluations in warm start
   %Prob.WarmStart = 1;
   %Result = rbfSolve(Prob);
   %PrintResult(Result);
end
if Solvers(3) | Solvers(4)
   Prob.PriLevOpt = 0;
   Prob.optParam.IterPrint = 1;     % Print one line each iteration
   %Prob.optParam.MaxFunc   = 400;   % Max function evaluations in ego

   if Solvers(4)
      Prob.CGO.globalSolver = 'glcCluster';
      %Prob.GO.maxFunc1 = 401;
      %Prob.GO.maxFunc2 = 200;
   else
      Prob.CGO.globalSolver = 'glcFast';
   end

   % Print EpsF
   EpsF=Prob.optParam.eps_f

%profile on
   Result = tomRun('ego',Prob,2);

%profile report
   %Prob.optParam.MaxFunc   = 20;  % Max function evaluations in warm start
   %Prob.WarmStart = 1;
   %Result = rbfSolve(Prob);
   %PrintResult(Result);
end

if Solvers(5)
   fprintf('\nRun Tomlab / glbFast\n\n');
   Prob.optParam.MaxFunc=40000;
   Prob.optParam.MaxIter=10000;
   Prob.optParam.IterPrint = 0;
   Result = tomRun('glbFast',Prob,2);
end
if Solvers(6)
   fprintf('\nRun Tomlab / glcFast\n\n');
   Prob.optParam.MaxFunc=40000;
   Prob.optParam.MaxIter=10000;
   Prob.optParam.IterPrint = 0;
   Result = tomRun('glcFast',Prob,2);
end

if Solvers(7)
   fprintf('\nRun Tomlab / glcCluster\n\n');
   Prob.optParam.MaxFunc=20000;
   Prob.optParam.MaxIter=1000;
   Prob.optParam.IterPrint = 0;
   Prob.PriLevOpt = 1;
   %Prob.GO.maxFunc1 = n*300+1;
   %Prob.GO.maxFunc2 = n*200;
   %Prob.GO.maxFunc1 = n*5000+1;
   %Prob.GO.maxFunc2 = n*200;
   %Prob.GO.localSolver = 'conSolve';
   %Prob.optParam.fGoal = -1E300;
   %Prob.MIP.fGoal = -1E300;
   %Prob.f_opt = -1E-300;

%profile on
   Result = tomRun('glcCluster',Prob,2);
%profile report

   %tic
   %%Result = glcCluster(Prob,11,10);
   %toc
   %fprintf('\n');
   % PrintResult also prints the CPU time
   %PrintResult(Result);
end

if Solvers(8)
   fprintf('\nRun Tomlab / glbSolve\n\n');
   Prob.PriLevOpt = 0;
   Prob.optParam.MaxFunc=20000;
   Prob.optParam.MaxFunc=5000;
   Prob.optParam.MaxIter=1000;
   Prob.optParam.MaxIter=500;
   Prob.optParam.IterPrint = 0;
   profile on
   tic
   Result = glbSolve(Prob);
   fprintf('\n');
   toc
   profile report
   % PrintResult also prints the CPU time
   PrintResult(Result);
end
if Solvers(9)
   fprintf('\nRun Tomlab / glcSolve\n\n');
   Prob.PriLevOpt = 0;
   Prob.optParam.MaxFunc=20000;
   Prob.optParam.MaxFunc=5000;
   Prob.optParam.MaxIter=1000;
   Prob.optParam.MaxIter=500;
   Prob.optParam.IterPrint = 0;
   tic
   Result = glcSolve(Prob);
   fprintf('\n');
   toc
   % PrintResult also prints the CPU time
   PrintResult(Result);
end

diary off

% MODIFICATION LOG:
%
% 080410 hkh Written 
