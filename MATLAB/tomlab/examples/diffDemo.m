% Demonstration examples for numerical differentiation
%
% The demonstration functions are
%
% diffHelp   Help on the differentiation examples
% diff1Demo  Solve Rosenbrocks banana (RB) with MAD automatic differentiation
% diff2Demo  Solve RB only
% diff3Demo  Spline example with Matlab spline
% diff4Demo  Illustrates how to set nonzeros in HessPattern

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Dec 27, 2001.   Last modified Apr 7, 2004.

function diffDemo

%----------------------------------------------
%This is the file diffDemo.m in tomlab\examples
%----------------------------------------------

diff1Demo
diff2Demo
diff3Demo
diff4Demo

% ------------------------------------------------------------------------
function diff1Demo
% ------------------------------------------------------------------------

format compact
fprintf('======================================================\n');

echo on

% Automatic Differentiation example}, if MAD TB is installed
Prob            = probInit('uc_prob',1);

madinitglobals;              % Initiate MAD variables

Prob.Solver.Alg = 2;         % Use the safeguarded standard BFGS
Prob.ADObj      = 1;         % Use Automatic Differentiation for objective

% Solve using tomRun driver routine
Result = tomRun('nlpSolve',Prob,2);

   
% ------------------------------------------------------------------------
function diff2Demo
% ------------------------------------------------------------------------

% Finite differentiation using the FD algorithm

fprintf('\n');
fprintf('------------------------------------------------------------\n');
disp('Solve Rosenbrocks banana problem with different Prob.GradTolg values')
disp('First use default values - Optimal value should be exactly 0')
fprintf('------------------------------------------------------------\n');
fprintf('\n');
Prob            = probInit('uc_prob',1);
Prob.Solver.Alg = 2;             
Prob.NumDiff    = 1;             % Use the fdng routine with the FD algorithm.
Result          = tomRun('ucSolve',Prob,1);
pause

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf('\n');
% Change the tolerances used by algorithm FD
disp('Change Prob.GradTolg to values 1E-5,1E-6')
disp('Results are a bit worse')
fprintf('------------------------------------------------------------\n');

Prob.GradTolg   = [1E-5; 1E-6];  % Change the predefined step size
Result          = tomRun('ucSolve',Prob,1);

% The change leads to worse accuracy

pause

fprintf('\n');
fprintf('\n');
fprintf('------------------------------------------------------------\n');
% Instead let an algorithm determine the best possible GradTolg
% It needs some extra f(x) evalutations, but the result is much better.

disp('Change Prob.GradTolg to -1, to let algorithm determine step length')
disp('Results improve');
fprintf('------------------------------------------------------------\n');
fprintf('\n');
Prob.GradTolg         = -1;  % A negative number demands that the step length
                             % algorithm is to be used at the initial point
% Avoid setting GradTolg empty, then instead Prob.optParam.DiffInt is used.

Result          = tomRun('ucSolve',Prob,1);

% Now the result is perfect, very close to the optimal == 0.

fprintf('------------------------------------------------------------\n');
fprintf('step size : %15.8e %15.8e\n',Result.Prob.GradTolg);
fprintf('------------------------------------------------------------\n');
pause

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf('\n');
Prob.NumDiff    = 5;             % Use the complex variable technique

fprintf('------------------------------------------------------------\n');
disp('Use Complex variable technique, NumDiff=5, gives perfect result!')
disp('The advantage is that it is insensitive to the choice of step length');
fprintf('------------------------------------------------------------\n');

Result          = tomRun('ucSolve',Prob,1);

% When it works, like in this case, it gives absolutely perfect result.

fprintf('\n');
fprintf('\n');
fprintf('------------------------------------------------------------\n');
disp('When it works, like in this case, it gives absolutely perfect result.')
fprintf('f_k with 30 decimals: %40.30f\n',Result.f_k);
fprintf('------------------------------------------------------------\n');

pause


function diff3Demo

   % Spline example with Matlab spline
   fprintf('------------------------------------------------------------\n');
   disp('Use the Matlab spline routine to estimate the least squares Jacobian')
   probFile        = 'ls_prob';
   fprintf('------------------------------------------------------------\n');
   P               = 1;
   Prob            = probInit(probFile, P);
   Prob.Solver.Alg = 0;             % Use the default algorithm in clsSolve
   Prob.NumDiff    = 2;             % Use the Matlab spline routine
   Result          = clsSolve(Prob);
   PrintResult(Result,2);
   pause
echo off
pause

function diff4Demo
% Finite differentiation using the FD algorithm
%
% This example illustrates how to set nonzeros in HessPattern
% to tell TOMLAB which elements to estimate in the Hessian.
% All elements in Hessian with corresponding zeros in HessPattern are
% set to 0, and no numerical estimate is needed.
%
% This saves very much time for large problems
% In a similar way, Prob.ConsPattern is set for the constraint gradient
% matrix for the nonlinear constraints, and Prob.JacPattern for the
% Jacobian matrix in nonlinear least squares problems.
disp('Solve constrained problem setting nonzeros in HessPattern for');
disp('Hessian of objective and estimating numerically using in');
disp('first case finite differences');
fprintf('\n');
pause

probFile                 = 'con_prob';
P                        = 12;
Prob                     = probInit(probFile, P);
Prob.Solver.Alg          = 1;             
Prob.HessPattern         = sparse([ones(6,6), zeros(6,6);zeros(6,12)]);

% Note that if setting Prob.NumDiff = 1, also the gradient would be
% estimated with numerical differences, which is not recommended.
% Estimating two levels of derivatives is an ill-conditioned process.
% Setting Prob.NumDiff = -1, only the Hessian is estimated

% Use qpSolve in base module for QP subproblems
Prob.SolverQP            = 'qpSolve';  

Prob.NumDiff             = -1;  % Use the fdng routine with the FD algorithm.
% Prob.optParam.IterPrint  = 1; % Displays one line per iteration
Result                   = tomRun('nlpSolve',Prob,1);
fprintf('\n');
fprintf('\n');
disp('Setting Prob.NumDiff = -11 makes Tomlab analyze HessPattern');
disp('and use a  faster algorithm for large-scale problems.');
disp('In this small-scale case it is no advantage,');
disp('it just takes more CPU-time');
pause
Prob.NumDiff             = -11;  % Use the fdng routine with the FD algorithm.
Result                   = tomRun('nlpSolve',Prob,1);
fprintf('\n');
fprintf('\n');

disp('Run the same problem estimating Hessian with Matlab Spline');
disp('Needs more gradient calculations because it is principle');
disp('smooth central differences to obtain the numerical Hessian.');
disp('If the problem is noisy, then this method is recommended');
fprintf('\n');
pause
Prob.NumDiff             = -2;  % Matlab spline
Result                   = tomRun('nlpSolve',Prob,1);
pause
