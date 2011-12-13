% demouc
% 
% Demonstration of the solution of unconstrained optimization problems (uc) 
% in TOMLAB.
%
% demouc;
%
% demouc makes echo on, displays comments about the solution
% process and pause statements
%
% If instead doing:
% PAUS=0;
% demouc
%
% Then demouc runs without pause statements and display of comments
% Suitable to check that the system is working properly

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Written June 30, 1999. Last modified Apr 10, 2004.

%eval('echo on');

echo('on')

%# call uc_prob

if ~exist('PAUS')
   PAUS=1;
end

if PAUS
   echo('on')
end

% Some different tests on the Rosenbrocks banana function
%

if PAUS, pause, end

% Run ucSolve with default values, i.e. BFGS quasi-Newton algorithm
% When only given uc_f, ucSolve will estimate the gradient with
% numerical differences

% First define the problem with a call to probAssign
%
Prob=probAssign('uc',[],[],'RB test',[-1.2 1]);

%
% Then tell TOMLAB which file computes the function value 
Prob=tomFiles(Prob,'uc_f');
%
% Now call the driver routine, tomRun.
% tomRun must know which solver to use (ucSolve) and the problem
% definition structure Prob
%

if PAUS, pause, end

tomRun('ucSolve',Prob,[],2);

if PAUS, pause, end

%
% If adding the gradient vector, the number of function value evaluations
% (FuncEv) are greatly reduced. 
%
% Also better accuracy is obtained, showed by the gradient being closer to 0 
% and the function value being smaller, closer to the optimal 0 value
%
% Call tomFiles once again with the name of the gradient routine

Prob=tomFiles(Prob,'uc_f','uc_g');

if PAUS, pause, end

tomRun('ucSolve',Prob,[],2);

if PAUS, pause, end

%
% It is also possible to run Newtons method, but the matrix of second
% derivatives, the Hessian matrix, must be estimated with numerical
% differences. 
%
% Set Prob.Solver.Alg==1 to get Newtons method.

Prob.Solver.Alg=1;

if PAUS, pause, end

tomRun('ucSolve',Prob,[],2);


% The number of function evaluations (FuncEv) and the number of iterations
% (Iter) are now lower, but the number of gradient evalutations (GradEv) are
% very much higher than with the BFGS method. 
%
% If also adding the Hessian matrix, the number of gradient value evaluations
% are greatly reduced. 
%
% Call tomFiles to also set the function to compute the Hessian. 

Prob=tomFiles(Prob,'uc_f','uc_g','uc_H');

% Now we can run Newtons method to get the fastest possible method for
% this type of problem

if PAUS, pause, end

tomRun('ucSolve',Prob,[],2);

if PAUS, pause, end

% We can also directly call a solver, without using tomRun.
% This is possible if we have been using probAssign/tomFiles to define
% the Prob struct.
%
% We can use some routine for general constrained problems as well.
%
% We try conSolve
%

Result=conSolve(Prob);

if PAUS, pause, end

% Note that now we do not get any output (only if the print level in
% the optimization routine is high (Prob.PriLevOpt > 5)
%
% We can get the normal output with a call to the routine PrintResult
% Normal print level is 2, we now try the maximum (3)
%

PrintResult(Result,3) 

if PAUS, pause, end


% SETUP IS USED
%
% When running many optimization problems it is best to use a
% setup routine, and define many problems in one file.
%
% Predefined unconstrained problems are available in uc_f,uc_g and uc_H.
% The Init File defining the test problems is uc_prob.
%
% To use the Init File approach, first define which file to be used,
% and the first problem to be solved

Prob=probInit('uc_prob',1);


% Then it possible to solve the problem using a driver routine or
% by a direct call to the solver
%
% First solve the problem using sTrustr. This routine is an alternative
% to ucSolve for unconstrained problems

if PAUS, pause, end

Result=tomRun('sTrustr',Prob,[],2);


if PAUS, pause, end

% Use Structured Trust Region with BFGS Update

Prob.Solver.Alg=1;


Result=tomRun('sTrustr',Prob,[],2);


% Numerical differences on Hessian, using complex variable method

Prob.NumDiff=-5;

Result=tomRun('sTrustr',Prob,[],2);

if PAUS, pause, end


% Numerical differences on both gradient and Hessian

Prob.NumDiff=1;

Result=tomRun('sTrustr',Prob,[],2);

if PAUS, pause, end


% Numerical differences on both gradient and Hessian with the complex 
% variable method

Prob.NumDiff=5;

Result=tomRun('sTrustr',Prob,[],2);

if PAUS, pause, end

% Numerical differences on Hessian
% Change of tolerances
% Prob.GradTolg (Prob.GradTolJ), and Prob.GradTolH

Prob.NumDiff=-1;
Prob.GradTolH=1E-10;
Prob.GradTolJ=1E-10;
% Used for Jacobian computations: Prob.GradTolJ=1E-10;

Result=tomRun('sTrustr',Prob,[],2);

if PAUS, pause, end


% Change the numerical tolerance to much higher value
% Also NumDiff=1 means that BOTH the gradient and the Hessian should be
% approximated

Prob.NumDiff=1;
Prob.GradTolH=1E-5;

Result=tomRun('sTrustr',Prob,[],2);

if PAUS, pause, end



% Try all six different search methods in ucSolve for Newtons method
%
% The method number 0-5 is set in Prob.optParam.method
% The first three, SVD, LU and LU with pivoting are recommended methods,
% The others are included for illustrative purposes, but might work OK.
% They all work well for this test example.
%
% Method four fools Matlab to use QR, the fifth method is using the
% built in inverse, and the sixth is using an explicitely computed inverse.

Prob=probInit('uc_prob',1);
for i=0:5
    Prob.Solver.Alg      = 1;
    Prob.optParam.method = i;
    Prob.PriLev          = 2;
    Result               = tomRun('ucSolve',Prob,[],2);
end
