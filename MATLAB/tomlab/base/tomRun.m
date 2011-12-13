% tomRun - General driver routine for TOMLAB
%
% If using the TOMLAB format (TQ), call with:
%
%   function [Result] = tomRun(Solver, Prob, PriLev, ask);
%
% NOTE! The call Result = tomRun(Solver, Prob, [], 2); will be interpreted
% as (for compatability reasons) PriLev = 2, ask = []
% (reverse order of parameters as in Tomlab 1.0 - 4.4)
%
% The following call will also work (similar 6-input format as below)
%
%   function [Result] = tomRun(Solver, [], [], Prob, PriLev, ask);
%
% If using the TOMLAB Init File format, call with:
%
%   function [Result] = tomRun(Solver, probFile, probNumber, Prob, PriLev, ask);
%
% A third alternative is the call
%
%   function [Result] = tomRun(Solver, probType, probNumber, Prob, PriLev, ask);
%
% Then the default file for problems of type probType is used.
%
% If calling with tomRun; (no arguments), a list of available solvers is given
%
% If calling with tomRun(probType);
%    a list of available solvers for probType is given
%
% tomRun checks if the second argument is a string, a structure, a number,
% or is empty, to determine which input format is used.
%
% if isempty(Solver),   tomRun is using the TOMLAB default solver for the
%                       problem type as default (calling GetSolver)
% if isempty(probFile), tomRun is using con_prob as default
%
%
% A problem available in the TOMLAB Init File format is defined using a call
%           Prob=probInit(probFile, probNumber, ask, Prob)
%
% INPUT: (if [] is given or less parameters are given, default values are used)
%
% Solver     The name of the solver that should be used to optimize the problem.
%            If the Solver may run several different optimization algorithms,
%            then the values of Prob.Solver.Alg and Prob.Solver.SubAlg
%            determines which algorithm.
%            The Solver name is put in Prob.Solver.Name
%
% probFile   User problem initialization file.
%
% probNumber Problem number in probFile.
%            If empty of left out, either probNumber=Prob.P (if set) or
%            otherwise probNumber=1.
%            When calling the probFile with probNumber=0, probFile must
%            return a string matrix with the names of the problems defined.
%
% Prob       Problem structure. Either define the structure with the
%            call:  Prob=probInit(probFile,probNumber,ask);
%            or set in Prob the parameters with special values.
%            See the manual for a description of the Prob structure
%
%            P        probNumber (YOU MUST SET THIS, IF SETTING PROB AS INPUT!)
%
%            Examples of other fields to set:
%            probFile probFile
%            uP       User problem parameters, which you can use when computing
%                     the functions. See the variable ask .
%            optParam A substructure with optimization parameters. Default value
%                     from optParamSet(Solver,probType)
%            x_0      Starting point for the problem.
% PriLev     Print level when displaying the result of the optimization in
%            routine PrintResult.
%            =0 No output, =1 Final result, shorter version,
%            =2 Final result, longer version, =3 All output
%            If isempty(PriLev), Prob.PriLev is used, if nonempty
%
%            The printing level in the optimization solver is controlled
%            by setting the parameter Prob.PriLevOpt
%
% ask        If ask>=1: ask questions in probFile. If ask = 0: use defaults
%            If ask < 0: use values in user params uP if defined or defaults.
%            If isempty(ask): If length(uP) > 0, ask=-1, else ask=0
%
% OUTPUT:
% Result Structure with optimization results. See manual for description
%        Result.Prob holds the input Prob structure

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., $Release: 7.4.0$
% Written Mar 2, 1998.    Last modified April 2, 2010.

function Result = tomRun(Solver, Prob, P3, P4, P5, P6)

global n_f n_g n_H    % Count of function, gradient, Hessian evaluations
global n_c            % Count of constraint, constr.grad, and 2nd der evals
global n_r n_J        % Count of residual, Jacobian and 2nd der evaluations
global GlobalLevel
GlobalLevel = [];     % Initialize to zero depth for recursion
global optType probType
probType = [];
optType = [];

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print
if isempty(MAX_x)
   MAX_x=20;
end
if isempty(MAX_c)
   MAX_c=20;
end
if isempty(MAX_r)
   MAX_r=30;
end

if nargin < 2
   if nargin < 1
      for i=1:10
          fprintf('\nSolvers for problem type %d\n',i)
          z=SolverList(i);
          disp(z);
      end
   else
      if ~ischar(Solver)
         fprintf('\nSolvers for problem type %d\n',Solver)
         z=SolverList(max(1,min(10,round(Solver))));
         disp(z);
      else
         for i=1:10
             fprintf('\nSolvers for problem type %d\n',i)
             z=SolverList(i);
             disp(z);
         end
      end
   end
   return
end

if nargin < 6
    P6=[];
    if nargin < 5
        P5=[];
        if nargin < 4
            P4=[];
            if nargin < 3
                P3=[];
            end
        end
    end
end

if isstruct(Prob)
   TQ         = 1;
   if nargin == 4 && isempty(P3)
      PriLev     = P4;
      ask        = [];
   else
      PriLev     = P3;
      ask        = P4;
   end
   optType    = Prob.probType;
   probFile   = 0;  % Signal that a Prob struct is used
   probNumber = 1;
elseif ischar(Prob)
   TQ         = 0;
   probFile   = Prob;
   probNumber = P3;
   Prob       = P4;
   PriLev     = P5;
   ask        = P6;
   if isstruct(Prob)
      if isfield(Prob,'probType')
         optType    = Prob.probType;
      else
         optType    = 3;
      end
   else
      Prob = [];
      optType    = 3;
   end
elseif isempty(Prob)
   % Assume 4th argument is the Prob structure
   if nargin < 4
      error('When 2nd argument empty, tomRun needs 4 arguments');
   end
   if ~isstruct(P4)
      error('When 2nd argument empty, 4th argument must be a Prob structure');
   end
   TQ         = 1;
   Prob       = P4;
   PriLev     = P5;
   ask        = P6;
   optType    = Prob.probType;
   probFile   = 0;  % Signal that a Prob struct is used
   probNumber = 1;
else
   TQ         = 0;
   probType   = Prob;
   [DataFile]=nameprob(probType,0);
   probFile   = DataFile(1,:);
   if isempty(P3)
      probNumber = 1;
   else
      probNumber = P3;
   end
   Prob       = P4;
   PriLev     = P5;
   ask        = P6;
   optType    = probType;
end

if TQ == 0
   if isstruct(Prob)
      if ~isempty(probNumber)
         Prob.P=probNumber;
      end
      x_0=Prob.x_0;  % If input x_0 has correct length, use this as start
      if isempty(probNumber)
         probNumber=Prob.P;
      end
      if isempty(probFile)
         probFile=Prob.probFile;
      else
         Prob.probFile=[];
      end
      Prob.CHECK=0;
   else
      x_0=[];
   end

   if isempty(probFile)
      [DataFile,NameFile,DefFile]=nameprob(optType);
      probFile=DataFile(DefFile,:);
   end

   if isempty(probNumber), probNumber=1; end
   
   if isempty(PriLev)
      if ~isempty(Prob.PriLev)
         PriLev=Prob.PriLev;
      else
          PriLev = 2;
      end
   end

   if isempty(ask)
      if isstruct(Prob)
         if isempty(Prob.uP)
            ask=0;
         else
            ask=-1;
         end
      else
         ask=0;
      end
   end

   if isempty(probFile), probFile='con_prob'; end
else
   x_0 = Prob.x_0;
   if isempty(PriLev)
      if ~isempty(Prob.PriLev)
         PriLev=Prob.PriLev;
      else
          PriLev = 0;
      end
   end
end

% Check if Prob structure already defined
MENU=Prob.MENU;
if Prob.P <= 0
   MENU=0;
end

if ~ischar(probFile)
   if probFile==0
      MENU=1; % Accept the Prob structure as it is
      probType = Prob.probType;
   end
end

if MENU==0 
   % Define problem in structure Prob
   probFile=deblank(probFile);
   Prob = probInit(probFile, probNumber(1), ask, Prob);
   probType = Prob.probType;
end

if isempty(probType)
   probType=optType;
end

Prob = mkbound(Prob);

if ~isempty(x_0)
   if length(x_0)==Prob.N && Prob.N ~=0
      Prob.x_0 = x_0(:); % Should be column vector. Use input starting values
   end
end

if isempty(Solver) 
   isType = checkType([],probType);
   Solver = GetSolver(isType,Prob.LargeScale);
else
   Solver=deblank(Solver);
end

Result   = [];
NLLSVars = 0;
GlobVars = 0;

[TomV,os,TV]=tomlabVersion;

switch lower(Solver)
 case 'nlpsolve' % SQP. Fletcher-Leyffer
   Result = nlpSolve(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=n_c;
   NLLSVars = 1;
   GlobVars = 1;
 case 'consolve' % Schittkowski Augmented Lagrangian SQP
   Result = conSolve(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=n_c;
   NLLSVars = 1;
   GlobVars = 1;
 case 'strustr' % Structured Trust Region
   % sTrustr runs a structured trust region algorithm 
   % (Conn/Gould/Sartenaer/Toint)
   Result = sTrustr(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=n_c;
 case 'clssolve' %General TOMLAB constrained LS, clsSolve:
   Result = clsSolve(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=0;
 case 'glbsolve'
   Result= glbSolve(Prob);
 case 'glbfast'
   Result= glbFast(Prob);
 case 'glbdirect'
   Result= glbDirectTL(Prob);
 case 'ego'
   Result= ego(Prob);
 case 'ego05'
   Result= ego05(Prob);
 case 'glcsolve'
   Result= glcSolve(Prob);
 case 'glcfast'
   Result= glcFast(Prob);
 case 'glcdirect'
  Result= glcDirectTL(Prob);
 case 'glccluster'
   Result= glcCluster(Prob);
 case 'minlpsolve' % General TOMLAB MINLP 
   Result = minlpSolve(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   NLLSVars = 1;
   GlobVars = 1;
 case 'filmint' % TOMLAB MINLP - FilMINT outer approximation algorithm
   Result = FilMINT(Prob);
 case 'rbfsolve'
   Result= rbfSolve(Prob);
 case 'arbfmip'
   Result= arbfmip(Prob);
 case 'arbf'
   Result= arbf(Prob);
 case 'dynrbf'
   Result= dynrbf(Prob);
 case 'strustr' % Structured Trust Region
   Result = sTrustr(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=n_c;
 case 'milpsolve' % TOMLAB general (MI)LP solver
   Prob = ProbCheck(Prob,'lpSimplex',7);
   Result = milpSolveTL(Prob);
 case 'multimin'
   Result = multiMin(Prob);
 case 'multiminlp'
   Result = multiMINLP(Prob);
 case 'sepminlp'
  Result = sepMINLP(Prob);
 case 'lpsimplex' % TOMLAB general LP solver
   Result = lpSimplex(Prob);
 case 'qpsolve' % TOMLAB general QP solver
   Result = qpSolve(Prob);
 case 'mipsolve' % TOMLAB general Branch and bound MIP solver
   Result = mipSolve(Prob);
 case 'cutplane' % TOMLAB cutting plane solver
   Result = cutplane(Prob);
 case 'ucsolve'     % TOMLAB unconstrained optimization solver
   Result = ucSolve(Prob);
   Result.FuncEv=n_f;
   Result.GradEv=n_g;
   Result.ConstrEv=0;
% -----------------------
% SOL Solvers
% -----------------------
 case 'lpopt' % Run LPOPT
   if TV(2)
      Prob   = ProbCheck(Prob,'lpopt',8);
      Result = lpoptTL(Prob);
   else
      fprintf('No valid license for the LPOPT solver\n');
   end
 case 'qpopt' % Run QPOPT
   if TV(2)
      Prob   = ProbCheck(Prob,'qpopt',2);
      Result = qpoptTL(Prob);
   else
      fprintf('No valid license for the QPOPT solver\n');
   end
 case {'sqopt','sqopt7'} % Run SQOPT
   if TV(4)
      Prob=ProbCheck(Prob,'sqopt',2);
      Result = sqoptTL(Prob);
   else
      fprintf('No valid license for the SQOPT solver\n');
   end
 case 'lssol' % Run LSSOL
   if TV(3)
      Prob=ProbCheck(Prob,'lssol',5);
      Result = lssolTL(Prob);
   else
      fprintf('No valid license for the LSSOL solver\n');
   end
 case 'nlssol' % Run NLSSOL
   if TV(3)
      Prob=ProbCheck(Prob,'nlssol',6);
      Result = nlssolTL(Prob);
   else
      fprintf('No valid license for the NLSSOL solver\n');
   end
 case 'minos'
   if TV(2)
      Prob   = ProbCheck(Prob,'minos',3);
      Result = minosTL(Prob);
   else
      fprintf('No valid license for the MINOS solver\n');
   end
   NLLSVars = 1;
 case 'lp-minos'
   if TV(2)
      Prob   = ProbCheck(Prob,'minos',8);
      Result = minoslpTL(Prob);
   else
      fprintf('No valid license for the MINOS solver\n');
   end
 case 'qp-minos'
   if TV(2)
      Prob   = ProbCheck(Prob,'minos',2);
      Result = minosqpTL(Prob);
   else
      fprintf('No valid license for the MINOS solver\n');
   end
 case 'npsol' % Run NPSOL
   if TV(3)
      Prob=ProbCheck(Prob,'npsol',3);
      Result = npsolTL(Prob);
   else
      fprintf('No valid license for the NPSOL solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;
 case {'snopt','snopt7','snopt6'}
    if TV(4)
       Prob=ProbCheck(Prob,'snopt',3);
       Result = snoptTL(Prob);
    else
      fprintf('No valid license for the SNOPT solver\n');
    end
    NLLSVars = 1;
    GlobVars = 1;
 case 'snopt8' % SNOPT 8
    if TV(4)
       Prob=ProbCheck(Prob,'snopt',3);
       Result = snopt8TL(Prob);
    else
      fprintf('No valid license for the SNOPT 8 solver\n');
    end
    NLLSVars = 1;
    GlobVars = 1;
% -----------------------
% Other Solvers
% -----------------------
 case 'tlsqr'      % Run mex interface to LSQR linear least squares
   Prob   = ProbCheck(Prob,'Tlsqr',5);
   Result = TlsqrTL(Prob);
 case 'ql'
   Prob   = ProbCheck(Prob,'ql',2);
   Result = qlTL(Prob);
 case 'qld'
   Prob   = ProbCheck(Prob,'qld',2);
   Result = qldTL(Prob);
 case 'lsei'      % Run mex interface to LSEI linear least squares
   Prob   = ProbCheck(Prob,'lsei',5);
   Result = lseiTL(Prob);
 case 'bqpd'
   if TV(7)
      Prob   = ProbCheck(Prob,'bqpd',2);
      Result = bqpdTL(Prob);
   else
      fprintf('No valid license for the BQPD solver\n');
   end
 case 'miqpbb'
   if TV(7)
      Prob   = ProbCheck(Prob,'miqpbb',11);
      Result = miqpBBTL(Prob);
   else
      fprintf('No valid license for the miqpBB solver\n');
   end
 case {'filtersqp', 'filsqp'}
   if TV(7)
      Prob   = ProbCheck(Prob,'filterSQP',3);
      Result = filterSQPTL(Prob);
   else
      fprintf('No valid license for the filterSQP solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;
 case 'minlpbb'
   if TV(7)
      Prob   = ProbCheck(Prob,'minlpBB',12);
      Result = minlpBBTL(Prob);
   else
      fprintf('No valid license for the minlpBB solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;
 case 'pensdp'
   if TV(6)
      Prob   = ProbCheck(Prob,'PENSDP',13);
      Result = pensdpTL(Prob);
   else
      fprintf('No valid license for the PENSDP solver\n');
   end
 case 'penbmi'
   if TV(10)
      Prob   = ProbCheck(Prob,'PENBMI',14);
      Result = penbmiTL(Prob);
   else
      fprintf('No valid license for the PENBMI solver\n');
   end
 case 'pennon'
   if TV(10)
      'TODO licensing PENNON'
      Prob   = ProbCheck(Prob,'PENBMI',14);
      Result = pennonTL(Prob);
   else
      fprintf('No valid license for the PENNON solver\n');
   end

   case {'nlpql', 'nlpqlp'}
   if TV(16) 
      Prob   = ProbCheck(Prob,'nlpql',3);
      Result = nlpqlTL(Prob);
   else
      fprintf('No valid license for the NLPQL solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;
 case {'dfnlp', 'dfnlpd'}
   if TV(16) 
      Prob   = ProbCheck(Prob,'dfnlp',3);
      Result = dfnlpTL(Prob);
   else
      fprintf('No valid license for the DFNLP solver\n');
   end
 case {'nlpjob'}
   if TV(16) 
      Prob   = ProbCheck(Prob,'nlpjob',3);
      Result = nlpjobTL(Prob);
   else
      fprintf('No valid license for the NLPJOB solver\n');
   end
 case {'misqp'}
    disp('NAG :: FIX LICENSING')
    Prob   = ProbCheck(Prob,'misqp',3);
    Result = misqpTL(Prob);
 case {'miql'}
    disp('NAG :: FIX LICENSING')
    Prob   = ProbCheck(Prob,'miql',11);
    Result = miqlTL(Prob);
 case 'pdco'
   Prob   = ProbCheck(Prob,'pdco',3);
   Result = pdcoTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'pdsco'
   Prob   = ProbCheck(Prob,'pdsco',3);
   Result = pdscoTL(Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'oqnlp'
   if TV(14)
      Prob    = ProbCheck(Prob,'oqnlp',12);
      Result  = oqnlpTL(Prob);
   else
      fprintf('No valid license for the OQNLP solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;
 case 'msnlp'
   if TV(22)
      Prob    = ProbCheck(Prob,'msnlp',12);
      Result  = msnlpTL(Prob);
   else
      fprintf('No valid license for the MSNLP solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;   
 case 'lsgrg2'
   if TV(22)
      Prob    = ProbCheck(Prob,'lsgrg2',3);
      Result  = lsgrg2TL(Prob);
   else
      fprintf('No valid license for the LSGRG2 solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;   
 case {'knitro'}
   if TV(11) 
      Prob    = ProbCheck(Prob,'knitro',3);
      Result  = knitroTL(Prob);
   else
      fprintf('No valid license for the KNITRO solver\n');
   end  
   NLLSVars = 1;
   GlobVars = 1;   
 case {'conopt'}
   if TV(12)
      Prob = ProbCheck(Prob,'conopt',3);
      Result = conoptTL(Prob);
   else
      fprintf('No valid license for the CONOPT solver\n');
   end
   NLLSVars = 1;
   GlobVars = 1;
% -----------------------
% Xpress-MP solver
% -----------------------
 case {'xpress-mp','xpress','xpressmp'} % Run Xpress-MP
   if TV(8)
      Prob=ProbCheck(Prob,'xpress-mp',11);
      Result = xpressTL(Prob);
   else
      fprintf('No valid license for the Xpress-MP solver\n');
   end
% -----------------------
% CPLEX Solver
% -----------------------
 case 'cplex' % Run CPLEX
   if TV(9)
      Prob=ProbCheck(Prob,'CPLEX',11);
      Result = cplexTL(Prob);
   else
      fprintf('No valid license for the CPLEX solver\n');
   end
 case 'cplex11', % Run CPLEX 11
   if TV(9)
      Prob=ProbCheck(Prob,'CPLEX',11);
      Result = cplex11TL(Prob);
   else
      fprintf('No valid license for the CPLEX solver\n');
   end
      
% -----------------------
% GUROBI Solver
% -----------------------
 case 'gurobi' % Run GUROBI
   if TV(32)
      Prob=ProbCheck(Prob,'GUROBI',7);
      Result = gurobiTL(Prob);
   else
      fprintf('No valid license for the GUROBI solver\n');
   end
% ------------------------
% XA MEX file interface
% ------------------------
 case 'xa'
   if TV(15)
      Prob = ProbCheck(Prob,'xa',2);
      Result = xaTL(Prob);
   else
      fprintf('No valid license for the XA solver\n');
   end
% -----------------------
% Opt tbx 1.x interfaces
% -----------------------
 case 'constr'
   Result=opt15Run('CONSTR',Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'fmins'
   Result=opt15Run('FMINS',Prob);
 case 'leastsq'
   Result=opt15Run('LEASTSQ',Prob);
 case 'lp'
   Result=opt15Run('LP',Prob);
 case 'qp'
   Result=opt15Run('QP',Prob);
 case 'fminu'
   Result=opt15Run('FMINU',Prob);
% -----------------------
% Opt tbx 2.x interfaces
% -----------------------
 case 'fminunc'
   Result=opt20Run('FMINUNC',Prob);
 case 'fmincon'
   Result=opt20Run('FMINCON',Prob);
   NLLSVars = 1;
   GlobVars = 1;
 case 'fminsearch'
   Result=opt20Run('FMINSEARCH',Prob);
 case 'lsqnonlin'
   Result=opt20Run('LSQNONLIN',Prob);
 case 'lsqlin'
   Result=opt20Run('LSQLIN',Prob);
 case 'linprog'
   Result=opt20Run('LINPROG',Prob);
 case 'quadprog'
   Result=opt20Run('QUADPROG',Prob);
 case {'gp', 'coplgp'}
   Prob=ProbCheck(Prob,'gp',checkType('gp'));
   Result=coplgpTL(Prob);
 case 'geno'
   Prob=ProbCheck(Prob,'geno',checkType('glc'));
   Result=genoTL(Prob); 
 case 'lgo'
   Prob=ProbCheck(Prob,'lgo',12);
   Result=lgoTL(Prob);
 case 'slssolve'
   Result=slsSolve(Prob,PriLev);
   PriLev = 0;
 case 'infsolve'
   Result=infSolve(Prob,PriLev);
   PriLev = 0;
 case 'inflinsolve'
   Result=infLinSolve(Prob,PriLev);
   PriLev = 0;
 case 'goalsolve'
   Result=goalSolve(Prob,PriLev);
   PriLev = 0;
 case 'l1solve'
   Result=L1Solve(Prob,PriLev);
   PriLev = 0;
 case 'l1linsolve'
   Result=L1LinSolve(Prob,PriLev);
   PriLev = 0;
 case 'linratsolve'
   Result=linRatSolve(Prob,PriLev);
   PriLev = 0;
 case 'expsolve'
   Result=expSolve(Prob,PriLev);
   PriLev = 0;
 otherwise
   fprintf('Solver %s',Solver);
   fprintf(' NOT found\n');
   error('Illegal solver algorithm!')
end

if isempty(Result), EmptyResult(Solver), end

if any(Prob.probType==[4 5 6])
   if ~isempty(Result)
      Result.ResEv=n_r;
      Result.JacEv=n_J;
   end
end

if isfield(Prob,'AMPL')
    if ~isempty(Prob.AMPL)
        Result = postSolveAMPL(Result);
    end
end

% HKH CHECK THIS CAREFULLY, only used in mex routines before
if NLLSVars
   if any(Prob.probType==[4 6 11])
      Result.ResEv=n_r;
      Result.JacEv=n_J;
   end
end
if GlobVars
   if any(Prob.probType==[1 2 3 10])
      Result.FuncEv=n_f;
      Result.GradEv=n_g;
      Result.HessEv=n_H;
      Result.ConstrEv=n_c;
   elseif any(Prob.probType==[1 9])
      Result.FuncEv=n_f;
      Result.GradEv=n_g;
      Result.HessEv=n_H;
   elseif any(Prob.probType==[4 6 11])
      Result.FuncEv=n_f;
      Result.GradEv=n_g;
      Result.HessEv=n_H;
      Result.ConstrEv=n_c;
   end
end

PrintResult(Result,PriLev);

function EmptyResult(Code)
fprintf('\n\n')
fprintf('tomRun: ');
fprintf('ERROR! %s solver returned empty Result',Code);
fprintf('\n\n')
fprintf('Run startup to set paths\n');
fprintf('\n\n')

% MODIFICATION LOG:
%
% 981013  hkh  Added call to iniSolve and endSolve. Now call MEX with mexRun
%              Deleted 2nd output parameter Prob. Instead in Result.Prob.
%              PrintResult now function, called with Result and Pri.
% 981022  hkh  Delete global p_f, p_g etc.
% 981026  hkh  Redefine input to ???Run files.
% 981102  hkh  Set an initial value on Prob.P. Matlab is giving warnings
%              when setting empty, so set the more dangerous value of -1.
% 981108  hkh  Changed PriLev levels from 0-3
% 981110  hkh  Improve comments, discussing default values
%              Use Prob.probFile, if set in Prob and empty as argument
%              Change to use globals MAX_x, MAX_c, MAX_r
% 981118  hkh  Add call to iniSolve to clear globals.
% 981129  hkh  Add global variable xGUI
%              Write messages to GUI window if xGUI true. 
% 990629  hkh  Made it a function
% 000726  hkh  Add snopt solver
% 000820  hkh  Use mexSOL for SOL solvers, change call to tomRun
% 000820  hkh  Update for v3.0
% 000925  hkh  Revision, remove optParam
% 001106  hkh  General use of tomRun, and no other driver routines
% 010715  hkh  Adding glbFast
% 010815  hkh  Adding glcFast
% 011111  hkh  Adding glcCluster and rbfSolve
% 011204  hkh  Adding lsqlin
% 020512  hkh  Adding Dundee QP solver bqpd
% 020621  hkh  Adding Dundee MIQP solver MIQPbb
% 020630  hkh  Adding Dundee solvers MINLPbb, filterSQP; and PENSDP
% 020702  hkh  Adding CPLEX
% 030110  hkh  Do preSolve if Prob.optParam.PreSolve
% 030116  hkh  Change lsqr to Tlsqr
% 030123  hkh  Add pdco and pdsco
% 030129  hkh  Remove empty varargin in call to opt20Run and opt15Run
% 030129  ango Edit names, Dundee solvers
% 030206  hkh  Change filSQP to filterSQP in mexRun call
% 030317  ango Add OQNLP
% 030427  hkh  Add KS solver nlpql
% 030514  ango Add KNITRO
% 030603  ango Add CONOPT
% 030613  ango Add dfnlp
% 030709  ango Add nlpjob
% 030728  medv Added postSolve for AMPL problems.
% 031023  ango BARNLP, SPRNLP added
% 031029  ango XA added
% 031103  ango Synonyms for xpress-mp added: xpress, xpressmp
% 040101  hkh  NPOPT deleted
% 040102  hkh  Avoid tests, do them in ProbCheck
% 040101  hkh  Make EmptyResult routine, incl. missing dll check
% 040109  ango Add LGO
% 040111  hkh  Removed calls to iniSolve and endSolve, done in solvers and TL
% 040120  ango SNOPT 7 added
% 040419  med  Added PATH
% 041023  hkh  Added minlpSolve
% 041123  hkh  Default PriLev 2 for probFile input, otherwise 0 
% 041123  hkh  Revision to change order of PriLev and ask
% 041216  med  OQNLP and MSNLP added
% 050124  ango SQOPT7 added
% 050214  frhe BOEING solvers added
% 050310  frhe glbDirect added (new version of glbFast)
% 050323  med  Added LPSOLVE and moved lpSimplex to lpSolve2
% 050330  ango lpSolve2 is now lpSimplex
% 050422  ango lpSolve changed to milpSolve
% 050422  hkh  Add solvers infSolve,goalSolve,slsSolve,L1Solve,L1LinSolve
% 050602  hkh  Add solver coplgp (also OK with gp)
% 050602  hkh  Use checkType for GP and PATH
% 050605  hkh  Add snopt8, and use snopt and snopt6 to call snoptTL (snopt6)
% 050613  hkh  Wrong TV code for nlpql solvers
% 050801  med  isstr replaced by ischar
% 060201  med  Interface solver printing fixed
% 060715  hkh  Add multiMin
% 060822  ango Add GENO
% 070517  hkh  Add sepMINLP
% 071008  hkh  Now rbfSolve, arbf, arbfmip and dynrbf for /CGO
% 080606  med  Removed LineParamDet
% 080606  med  Cleaned the code
% 081112  ango sqopt7 replaces sqopt6
% 090228  med  socs removed
% 090508  med  gurobi added
% 090910  hkh  multiMINLP added
% 091118  hkh  FilMINT added
% 100402  hkh  miql added (also misqp)
