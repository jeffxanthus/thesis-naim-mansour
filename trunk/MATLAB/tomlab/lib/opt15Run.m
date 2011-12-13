% opt15Run: Code to run optimization solvers in MathWorks Optimization TB 1.5
%
% function Result = opt15Run(Solver, Prob)
%
% Solver is currently one of:
%    CONSTR
%    FMINS
%    FMINU
%    LP
%    QP
%    LEASTSQ
%
% The test on the name in Solver is not case-sensitive
%
% Prob   is the TOMLAB problem structure
%
% Result is the TOMLAB result structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.2.0$
% Written July 5, 1999.  Last modified Jun 6, 2008.

function Result = opt15Run(Solver, Prob)

global n_f n_r

if strcmpi(Solver,'FMINS') 
   Prob=ProbCheck(Prob,'FMINS',1);

   PriLev=Prob.PriLevOpt;

   Prob=iniSolve(Prob,checkType('uc'),0,0);

   Result=ResultDef(Prob);
   Result.Solver='fmins';
   Result.SolverAlgorithm='Nelder-Mead simplex algorithm in Optimization TB';

   optp=optimDef(Prob.optParam,Prob.LineParam,Prob);

   % MATLAB fmins simplex method

   if PriLev >= 0
      disp('Call MATLAB simplex routine fmins');
   end
   if optp(1)==1, optp(1)=0; end % optp(1)==1 displays each iter in fmins
   optp(5)=0;
   
  [x_k, optp]=fminsearch('nlp_f',Prob.x_0,optp,[],Prob);
      
  Result=endSolve(Prob,Result);
  Result.x_k=x_k;
  Result.f_k = optp(8);
  Result.FuncEv=optp(10);
  Result.f_0 = nlp_f(Result.Prob.x_0,Result.Prob);
  n_f = n_f - 1;  % Do not count this evalutation
  % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
  Result.xState=(x_k==Prob.x_L)+2*(x_k==Prob.x_U);
  Result.Inform=0;
  Result.ExitFlag=0;
else
  Result = opt1xRun(Solver, Prob);
end

if ~isempty(Result)
   Result=endSolve(Prob,Result);
   if isempty(Result.f_0)
      % Compute (x_0,f_0)
      nf=n_f;
      nr=n_r;
      if ~isempty(Result.x_0)
         Result.f_0 = nlp_f(Result.x_0,Result.Prob);
      end
      % Do not count this extra evaluation
      n_f=nf;
      n_r=nr;
   end
end

% MODIFICATION LOG:
%
% 981028  hkh  Add computation of (x_0,f_0)
% 981111  hkh  Use Result.x_0, not Result.Prob.x_0
% 981116  hkh  Change logic to count extra evaluation
% 990222  hkh  Only compute Result.f_0 if empty (not already computed)
% 990502  hkh  Add option to run pure LP in minos interface
% 030129  hkh  Remove varargin, always empty
% 040415  hkh  Make a correct call to iniSolve