% opt20Run: Code to run optimization solvers in MathWorks Optimization TB 2.0
%
% function Result = opt20Run(Solver, Prob)
%
% Solver is currently one of:
%    LINPROG
%    QUADPROG
%    FMINCON
%    FMINUNC
%    LSQNONLIN
%    FMINSEARCH
%    LSQLIN
%
% The test on the name in Solver is not case-sensitive
%
% Prob   is the TOMLAB problem structure
%
% Result is the TOMLAB result structure

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written July 5, 1999.  Last modified Aug 4, 2009.

function Result = opt20Run(Solver, Prob)

global n_f n_r n_g
global otxProb

PriLev=Prob.PriLevOpt;

if isfield(Prob,'Solver') 
   if isfield(Prob.Solver,'Alg') 
      alg=Prob.Solver.Alg;
   else
      alg=0;
   end
else
   alg=0;
end
if isempty(alg), alg=0; end
Prob.Solver.Alg=alg;

optstruct = optimset;

if strcmpi(Solver,'linprog') | strcmpi(Solver,'quadprog')
   % Define LP matrices and vectors from Prob structure
   if strcmpi(Solver(1),'l')
      Prob=ProbCheck(Prob,'LINPROG',8);
   else
      Prob=ProbCheck(Prob,'QUADPROG',2);
   end

   c=Prob.QP.c(:);

   A=Prob.A;
   b_L=Prob.b_L;
   b_U=Prob.b_U;
   x_L=Prob.x_L;
   x_U=Prob.x_U;
   n=max([length(c);length(x_L);length(x_U);size(A,2)]);
   if isempty(c)
      c=zeros(n,1);
   end

   bEqual=b_L==b_U;
   if strcmpi(Solver,'quadprog')
      F=Prob.QP.F;
   end
end


if strcmpi(Solver,'LINPROG')
   Prob=iniSolve(Prob,checkType('lp'),0,0);
   Result=ResultDef(Prob);
   Result.Solver='linprog';

   % Call MATLAB OPTIM:s new LP-routine

   
   Result.SolverAlgorithm='linprog';

   if PriLev > 0
      fprintf('\n\nCall MATLAB Optimization TB routine linprog\n');
   end

   if exist('linprog','file')
      if any(isnan(Prob.x_0))
         Prob.x_0=[];
      end
      checkx0(Prob.x_0,Prob.x_L,Prob.x_U);
   
      [AA,bb,meq] = cpTransf(Prob,2,0);
      [mm nn]=size(AA);
      cc=[c;zeros(nn-n,1)];
      if nn > n & ~isempty(Prob.x_0)
         Prob.x_0 = [Prob.x_0;bb(mm-(nn-n)+1:mm)];
      end
      if Prob.Solver.Alg==0
         X0=Prob.x_0;
         optstruct.LargeScale= 'off';
      else
         X0=[];
         optstruct.LargeScale= 'on';
      end
      if PriLev <= 0
         optstruct.Display= 'off';
      elseif PriLev == 0
         optstruct.Display= 'final';
      else
         optstruct.Display= 'iter';
      end
      optstruct.TolFun    = Prob.optParam.eps_f;
      optstruct.TolX      = Prob.optParam.eps_x;
      optstruct.MaxIter   = Prob.optParam.MaxIter;
      if meq > 0
         if xnargin('linprog') == 10
            [x_k, Result.f_k,ExitFlag,out,lam,R]=linprog(...
            cc,AA(meq+1:mm,:),bb(meq+1:mm),AA(1:meq,:),bb(1:meq),...
            Prob.x_L, Prob.x_U, X0, optstruct,Prob);
            Result.ExitText = R.ExitText;
         else
            [x_k, Result.f_k,ExitFlag,out,lam]=linprog(...
            cc,AA(meq+1:mm,:),bb(meq+1:mm),AA(1:meq,:),bb(1:meq),...
            Prob.x_L, Prob.x_U, X0, optstruct);
         end
      else
         if xnargin('linprog') == 10
            [x_k, Result.f_k,ExitFlag,out,lam,R]=linprog(...
            cc,AA,bb,[],[], Prob.x_L, Prob.x_U, X0, optstruct,Prob);
            Result.ExitText = R.ExitText;
         else
            [x_k, Result.f_k,ExitFlag,out,lam]=linprog(...
            cc,AA,bb,[],[], Prob.x_L, Prob.x_U, X0, optstruct);
         end
      end

      Result.Inform=ExitFlag;
      ExitFlag=max(-1,ExitFlag);
      switch ExitFlag
         case 0
           ExitFlag=1;  % Too many iterations
         case 1
           ExitFlag=0;  % OK
         case 2
           ExitFlag=0;  % OK
         case -1
           ExitFlag=3;  % Failure, but what type???
         case 10000
           ExitFlag=2;  % Unbounded ???
         case 10001
           ExitFlag=4;  % Infeasible ???
         case 10002
           ExitFlag=3;  % Rank problem ???
         case 10003
           ExitFlag=10; % Input errors ???
         otherwise
           ExitFlag=0;  % OK
      end

      Result.ExitFlag=ExitFlag;
      Result.v_k=[lam.eqlin;lam.ineqlin];
      x_k=x_k(1:n);
      Result.x_k=x_k;
      Result.x_0 = Prob.x_0;
      if ~isempty(Prob.x_0)
         Result.f_0 = c'*Prob.x_0;
      end

      xTol=Prob.optParam.xTol;
      B = double((x_k > x_L + xTol) & (x_k < x_U - xTol) & (x_L~=x_U));
      B(find(x_k >= x_U-xTol))=-1;
      Result.QP.B=B;

      Result.FuncEv=out.iterations;
      Result.ConstrEv=Result.FuncEv;
      Result.Iter=out.iterations;
      Result.SolverAlgorithm=out.algorithm;

      % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
      Result.xState=(x_k==Prob.x_L)+2*(x_k==Prob.x_U);
      Result=endSolve(Prob,Result);
   else
      fprintf('Can not find linprog.m.\n')
      fprintf('Have you got a license for Optimization Toolbox?\n')
      Result=[];
   end
elseif strcmpi(Solver,'QUADPROG')
   Prob=iniSolve(Prob,checkType('qp'),0,0);
   Result=ResultDef(Prob);
   Result.Solver='quadprog';

   % Call MATLAB OPTIM:s new QP-routine

   
   Result.SolverAlgorithm='quadprog';

   if PriLev > 0
      fprintf('\n\nCall MATLAB Optimization TB routine quadprog\n');
   end

   if exist('quadprog','file')
      if any(isnan(Prob.x_0))
         Prob.x_0=[];
      end
      checkx0(Prob.x_0,Prob.x_L,Prob.x_U);
   
      [AA,bb,meq] = cpTransf(Prob,2,0);
      [mm nn]=size(AA);
      cc=[c;zeros(nn-n,1)];
      if nn > n & ~isempty(Prob.x_0)
         Prob.x_0 = [Prob.x_0;bb(mm-(nn-n)+1:mm)];
      end
      if Prob.Solver.Alg==0
         X0=Prob.x_0;
         optstruct.LargeScale= 'off';
      else
         X0=[];
         optstruct.LargeScale= 'on';
      end
      if PriLev <= 0
         optstruct.Display= 'off';
      elseif PriLev == 0
         optstruct.Display= 'final';
      else
         optstruct.Display= 'iter';
      end
      optstruct.TolFun    = Prob.optParam.eps_f;
      optstruct.TolX      = Prob.optParam.eps_x;
      optstruct.MaxIter   = Prob.optParam.MaxIter;

      if xnargin('quadprog') == 11
         if meq > 0
            [x_k, Result.f_k,ExitFlag,out,lam]=quadprog(...
            F,cc,AA(meq+1:mm,:),bb(meq+1:mm),AA(1:meq,:),bb(1:meq),...
            Prob.x_L, Prob.x_U, X0, optstruct,Prob);
         else
            [x_k, Result.f_k,ExitFlag,out,lam]=quadprog(...
            F,cc,AA,bb,[],[], Prob.x_L, Prob.x_U, X0, optstruct,Prob);
         end
      else
         if meq > 0
            [x_k, Result.f_k,ExitFlag,out,lam]=quadprog(...
            F,cc,AA(meq+1:mm,:),bb(meq+1:mm),AA(1:meq,:),bb(1:meq),...
            Prob.x_L, Prob.x_U, X0, optstruct);
         else
            [x_k, Result.f_k,ExitFlag,out,lam]=quadprog(...
            F,cc,AA,bb,[],[], Prob.x_L, Prob.x_U, X0, optstruct);
         end
      end

      Result.Inform=ExitFlag;
      ExitFlag=max(-1,ExitFlag);
      switch ExitFlag
         case 0
           ExitFlag=1;  % Too many iterations
         case 1
           ExitFlag=0;  % OK
         case 2
           ExitFlag=0;  % OK
         case -1
           ExitFlag=3;  % Failure, but what type???
         case 10000
           ExitFlag=2;  % Unbounded ???
         case 10001
           ExitFlag=4;  % Infeasible ???
         case 10002
           ExitFlag=3;  % Rank problem ???
         case 10003
           ExitFlag=10; % Input errors ???
         otherwise
           ExitFlag=0;  % OK
      end

      Result.ExitFlag=ExitFlag;
      Result.v_k=[lam.eqlin;lam.ineqlin];
      x_k=x_k(1:n);
      Result.x_k=x_k;
      Result.x_0 = Prob.x_0;
      if ~isempty(Prob.x_0)
         Result.f_0 = qp_f(Prob.x_0,Prob);
      end
      xTol=Prob.optParam.xTol;
      B = (x_k > x_L + xTol) & (x_k < x_U - xTol) & (x_L~=x_U);
      B(find(x_k >= x_U-xTol))=-1;
      Result.QP.B=B;

      Result.FuncEv=out.iterations;
      Result.ConstrEv=Result.FuncEv;
      Result.Iter=out.iterations;
      Result.SolverAlgorithm=out.algorithm;

      % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
      Result.xState=(x_k==Prob.x_L)+2*(x_k==Prob.x_U);
      Result=endSolve(Prob,Result);
   else
      fprintf('Can not find quadprog.m.\n')
      fprintf('Have you got a license for Optimization Toolbox?\n')
      Result=[];
   end
elseif strcmpi(Solver,'FMINCON')

   % Call MATLAB OPTIM:s new con-routine

   Prob=ProbCheck(Prob,'FMINCON',3);

   Prob=iniSolve(Prob,checkType('con'),2,1);
   Result=ResultDef(Prob);
   Result.Solver='fmincon';
   Result.SolverAlgorithm='SQP algorithm in Optimization TB';

   if PriLev > 0
      fprintf('\n\nCall MATLAB Optimization TB routine fmincon\n');
   end

   if any(isnan(Prob.x_0))
      Prob.x_0=[];
   end
   checkx0(Prob.x_0,Prob.x_L,Prob.x_U);

   %[mA, Ax, bEqual, b_L, b_U, A] = LinConstr(Prob)

   %meq=sum(bEqual);

   [AA,bb,meq] = cpTransf(Prob,2,0);
   [mm nn]=size(AA);
   n=length(Prob.x_0);
   if nn > n & ~isempty(Prob.x_0)
      Prob.x_0 = [Prob.x_0;bb(mm-(nn-n)+1:mm)];
   end

   if Prob.Solver.Alg==0
      X0=Prob.x_0;
      optstruct.LargeScale= 'off';
   else
      X0=Prob.x_0;
      optstruct.LargeScale= 'on';
   end
   if PriLev <= 0
      optstruct.Display= 'off';
   elseif PriLev == 0
      optstruct.Display= 'final';
   else
      optstruct.Display= 'iter';
   end
   optstruct.TolFun      = Prob.optParam.eps_f;
   optstruct.TolX        = Prob.optParam.eps_x;
   optstruct.MaxIter     = Prob.optParam.MaxIter;
   optstruct.TolCon      = Prob.optParam.cTol;
   optstruct.MaxFunEvals = max(10000,Prob.optParam.MaxFunc);
   if isempty(Prob.FUNCS.g)
      optstruct.GradObj     = 'off';
   else
      optstruct.GradObj     = 'on';
   end
   if isempty(Prob.FUNCS.H)
      optstruct.Hessian     = 'off';
   else
      optstruct.Hessian     = 'on';
   end
   %optstruct.GradObj     = 'on';
   %optstruct.Hessian     = 'on';

   % NOT: Send info about functions to optim_fgH
   % Send info about functions to nlp_fgH
   Prob.FUNCSX.f   = Prob.FUNCS.f;
   Prob.FUNCSX.g   = Prob.FUNCS.g;
   Prob.FUNCSX.H   = Prob.FUNCS.H;
   % NOT: Send info about functions to optim_cdc
   % Send info about functions to nlp_cdceq
   Prob.FUNCSX.c   = Prob.FUNCS.c;
   Prob.FUNCSX.dc  = Prob.FUNCS.dc;
   Prob.FUNCSX.d2c = Prob.FUNCS.d2c;

   nA=nargout('fmincon');
   if nA > 7 
      % Call TOMLAB version of fmincon
      otxProb = Prob;
      Prob = rmfield(Prob,'TOMLAB');
      if meq > 0
         [x_k, Result.f_k, ExitFlag, out, lam, Result.g_k, Result.H_k,R]...
               = fmincon('nlp2_fgH',X0,...
         AA(meq+1:mm,:),bb(meq+1:mm),AA(1:meq,:),bb(1:meq),...
         Prob.x_L, Prob.x_U, 'nlp2_cdceq', optstruct, Prob);
      else
         [x_k, Result.f_k, ExitFlag, out, lam, Result.g_k, Result.H_k,R]...
               = fmincon('nlp2_fgH',X0, AA,bb,[],[],...
         Prob.x_L, Prob.x_U, 'nlp2_cdceq', optstruct, Prob);
      end
      Result.ExitText = R.ExitText;
      otxProb = [];
   else
      if meq > 0
         [x_k, Result.f_k, ExitFlag, out, lam, Result.g_k, Result.H_k]...
               = fmincon('nlp_fgH',X0,...
         AA(meq+1:mm,:),bb(meq+1:mm),AA(1:meq,:),bb(1:meq),...
         Prob.x_L, Prob.x_U, 'nlp_cdceq', optstruct, Prob);
      else
         [x_k, Result.f_k, ExitFlag, out, lam, Result.g_k, Result.H_k]...
               = fmincon('nlp_fgH',X0, AA,bb,[],[],...
         Prob.x_L, Prob.x_U, 'nlp_cdceq', optstruct, Prob);
      end
   end
   Result.Inform=ExitFlag;
   ExitFlag=max(-1,ExitFlag);
   switch ExitFlag
      case 0
        ExitFlag=1;  % Too many iterations
      case 1
        ExitFlag=0;  % OK
      case 2
        ExitFlag=0;  % OK
      case -1
        ExitFlag=3;  % Failure, but what type???
      case 10000
        ExitFlag=2;  % Unbounded ???
      case 10001
        ExitFlag=4;  % Infeasible ???
      case 10002
        ExitFlag=3;  % Rank problem ???
      case 10003
        ExitFlag=10; % Input errors ???
      otherwise
        ExitFlag=0;  % OK
   end

   Result.ExitFlag=ExitFlag;

   Result.v_k=[lam.lower;lam.upper;lam.eqlin;lam.ineqlin;...
               lam.eqnonlin;lam.ineqnonlin];
   Result.x_k=x_k;
   Result.x_0 = Prob.x_0;
   Result.FuncEv=out.funcCount;
   Result.ConstrEv=Result.FuncEv;
   Result.GradEv=n_g;
   Result.Iter=out.iterations;
   Result.SolverAlgorithm=out.algorithm;
   if ExitFlag <= 1
      % Must call function value routines to compute output that is missing
      Result.c_k  = nlp_c(x_k, Prob);
      Result.cJac = nlp_dc(x_k, Prob);
   end

   % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
   n = Prob.N;
   Result.xState=(x_k(1:n)==Prob.x_L)+2*(x_k(1:n)==Prob.x_U);
   Result=endSolve(Prob,Result);
   
elseif strcmpi(Solver,'FMINUNC')

   % Call MATLAB OPTIM:s new uc-routine

   Prob=ProbCheck(Prob,'FMINUNC',1);

   Prob=iniSolve(Prob,checkType('uc'),2,0);
   Result=ResultDef(Prob);
   Result.Solver='fminunc';

   if PriLev > 0
      fprintf('\n\nCall MATLAB Optimization TB 2.0 routine fminunc\n');
   end

   if any(isnan(Prob.x_0))
      Prob.x_0=[];
   end
   checkx0(Prob.x_0,Prob.x_L,Prob.x_U);

   if Prob.Solver.Alg==0
      X0=Prob.x_0;
      optstruct.LargeScale= 'off';
   else
      X0=Prob.x_0;
      optstruct.LargeScale= 'on';
   end
   if PriLev <= 0
      optstruct.Display= 'off';
   elseif PriLev == 0
      optstruct.Display= 'final';
   else
      optstruct.Display= 'iter';
   end
   optstruct.TolFun      = Prob.optParam.eps_f;
   optstruct.TolX        = Prob.optParam.eps_x;
   optstruct.MaxIter     = Prob.optParam.MaxIter;
   %optstruct.MaxFunEvals = 10000;
   %optstruct.MaxFunEvals = Prob.optParam.MaxFunc;
   optstruct.MaxFunEvals = max(10000,Prob.optParam.MaxFunc);
   if isempty(Prob.FUNCS.g)
      optstruct.GradObj     = 'off';
   else
      optstruct.GradObj     = 'on';
   end
   if isempty(Prob.FUNCS.H)
      optstruct.Hessian     = 'off';
   else
      optstruct.Hessian     = 'on';
   end

   % Send info about functions to optim_fgH
   Prob.FUNCSX.f   = Prob.FUNCS.f;
   Prob.FUNCSX.g   = Prob.FUNCS.g;
   Prob.FUNCSX.H   = Prob.FUNCS.H;


   if nargout('fminunc') == 7
      % Call TOMLAB version of fmincon
      otxProb = Prob;
      Prob = rmfield(Prob,'TOMLAB');
      [x_k, Result.f_k, ExitFlag, out, Result.g_k, Result.H_k,R]...
            = fminunc('nlp2_fgH',X0, optstruct, Prob);
      Result.ExitText = R.ExitText;
      otxProb = [];
   else
      [x_k, Result.f_k, ExitFlag, out, Result.g_k, Result.H_k]...
            = fminunc('nlp2_fgH',X0, optstruct, Prob);
   end

   Result.Inform=ExitFlag;
   ExitFlag=max(-1,ExitFlag);
   switch ExitFlag
      case 0
        ExitFlag=1;  % Too many iterations
      case 1
        ExitFlag=0;  % OK
      case 2
        ExitFlag=0;  % OK
      case -1
        ExitFlag=3;  % Failure, but what type???
      case 10000
        ExitFlag=2;  % Unbounded ???
      case 10001
        ExitFlag=4;  % Infeasible ???
      case 10002
        ExitFlag=3;  % Rank problem ???
      case 10003
        ExitFlag=10; % Input errors ???
      otherwise
        ExitFlag=0;  % OK
   end

   Result.ExitFlag=ExitFlag;

   Result.x_k=x_k;
   Result.x_0 = Prob.x_0;
   Result.FuncEv=out.funcCount;
   Result.GradEv=n_g;
   Result.Iter=out.iterations;
   Result.SolverAlgorithm=out.algorithm;

   % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
   Result.xState=(x_k==Prob.x_L)+2*(x_k==Prob.x_U);
   Result=endSolve(Prob,Result);

elseif strcmpi(Solver,'LSQNONLIN') 

   % Call MATLAB OPTIM:s new NLLS-routine

   Prob=ProbCheck(Prob,'LSQNONLIN',4);

   Prob=iniSolve(Prob,checkType('ls'),1,0);
   Result=ResultDef(Prob);
   Result.Solver='lsqnonlin';
   Result.SolverAlgorithm='NLLS algorithm in Optimization TB 2.0';

   if PriLev > 0
      fprintf('\n\nCall MATLAB Optimization TB routine lsqnonlin\n');
   end

   if any(isnan(Prob.x_0))
      Prob.x_0=[];
   end
   checkx0(Prob.x_0,Prob.x_L,Prob.x_U);
   if Prob.LineParam.LineAlg==1
      LineAlg='cubicpoly';
   else
      LineAlg='quadcubic';
   end
   if Prob.Solver.Alg==0
      % Levenberg-Marquardt
      alg='on';
      LargeScale='off';
   elseif Prob.Solver.Alg==1
      % Gauss Newton
      alg='off';
      LargeScale='off';
   else
      % Sparse
      alg='on';
      LargeScale='on';
   end
   if Prob.Solver.Method==0
      HessUpdate='bfgs';
   elseif Prob.Solver.Method==1
      HessUpdate='dfp';
   elseif Prob.Solver.Method==2
      HessUpdate='gillmurray';
   elseif Prob.Solver.Method==3
      HessUpdate='steepdesc';
   end
   if Prob.Solver.Alg==0
      optstruct.HessUpdate=HessUpdate;
      optstruct.LevenbergMarquardt=alg;
      optstruct.DiffMinChange=0.01*Prob.optParam.DiffInt;
      optstruct.DiffMaxChange=0.01;
   else
      optstruct.LargeScale= 'on';
   end

   if PriLev <= 0
      optstruct.Display= 'off';
   elseif PriLev == 0
      optstruct.Display= 'final';
   else
      optstruct.Display= 'iter';
   end
   optstruct.LargeScale     = LargeScale;
   optstruct.TolFun         = Prob.optParam.eps_f;
   optstruct.TolX           = Prob.optParam.eps_x;
   optstruct.MaxIter        = Prob.optParam.MaxIter;
   optstruct.TolCon         = Prob.optParam.cTol;
   optstruct.MaxFunEvals    = max(10000,Prob.optParam.MaxFunc);
   if isempty(Prob.FUNCS.J)
      optstruct.Jacobian    = 'off';
   else
      optstruct.Jacobian    = 'on';
   end
   optstruct.LineSearchType = LineAlg;

   % Send info about functions to ls2_rJ
   Prob.FUNCSX.r   = Prob.FUNCS.r;
   Prob.FUNCSX.J   = Prob.FUNCS.J;
   Prob.FUNCSX.d2r = Prob.FUNCS.d2r;
   Prob.FUNCSX.rJ  = Prob.FUNCS.rJ;
      
   nA=nargout('lsqnonlin');
   if nA > 7
      % Call TOMLAB version of lsqnonlin
      otxProb = Prob;
      Prob = rmfield(Prob,'TOMLAB');
      if strcmpi(LargeScale,'on')
         [x_k, Result.f_k, Result.r_k, ExitFlag, out, lam, ...
               Result.J_k, Result] = lsqnonlin('ls2_rJS',Prob.x_0, ...
               Prob.x_L, Prob.x_U, optstruct, Prob);
         Prob=Result.Prob;
         Prob.LargeScale=1;
      else
         % OPTIM TB does not handle bounds, but clsSolve does so we use it
         [x_k, Result.f_k, Result.r_k, ExitFlag, out, lam, ...
               Result.J_k, Result] = lsqnonlin('ls2_rJS',Prob.x_0, ...
               Prob.x_L, Prob.x_U, optstruct, Prob);
         Prob=Result.Prob;
         Prob.LargeScale=0;
      end
      otxProb = [];
   else
      if strcmpi(LargeScale,'on')
         % Large Scale algorithm can not handle fixed variables, crashes.
         eq=find(Prob.x_L==Prob.x_U);
         if ~isempty(eq)
            Prob.x_U(eq)=Prob.x_U(eq)+1E-10;
         end
         [x_k, Result.f_k, Result.r_k, ExitFlag, out, lam, ...
               Result.J_k] = lsqnonlin('ls2_rJS',Prob.x_0, ...
               Prob.x_L, Prob.x_U, optstruct, Prob);
         Prob.LargeScale=1;
      else
         % Does not handle bounds
         [x_k, Result.f_k, Result.r_k, ExitFlag, out, lam, Result.J_k]...
             = lsqnonlin('ls2_rJ',Prob.x_0, [], [], optstruct, ...
                          Prob);
         Prob.LargeScale=0;
      end
   end
   Result.Inform=ExitFlag;
   ExitFlag=max(-1,ExitFlag);
   switch ExitFlag
      case 0
        ExitFlag=1;  % Too many iterations
      case 1
        ExitFlag=0;  % OK
      case 2
        ExitFlag=0;  % OK
      case -1
        ExitFlag=3;  % Failure, but what type???
      case 10000
        ExitFlag=2;  % Unbounded ???
      case 10001
        ExitFlag=4;  % Infeasible ???
      case 10002
        ExitFlag=3;  % Rank problem ???
      case 10003
        ExitFlag=10; % Input errors ???
      otherwise
        ExitFlag=0;  % OK
   end

   Result.ExitFlag=ExitFlag;

   Result.f_k=0.5*Result.f_k; % OPT TB 2.0 returns r'*r, without 0.5 times.
   n = length(x_k);
   Result.v_k=[lam.lower;lam.upper];
   Result.x_k=x_k;
   Result.x_0 = Prob.x_0;
   Result.FuncEv=out.funcCount;
   Result.GradEv=n_g;
   Result.Iter=out.iterations;
   Result.g_k=Result.J_k'*Result.r_k;
   Result.SolverAlgorithm=out.algorithm;

   % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
   Result.xState=(x_k(1:n)==Prob.x_L)+2*(x_k(1:n)==Prob.x_U);
   Result=endSolve(Prob,Result);
   
elseif strcmpi(Solver,'FMINSEARCH')
   Prob=ProbCheck(Prob,'FMINSEARCH',1);

   Prob=iniSolve(Prob,checkType('uc'),0,0);
   Result=ResultDef(Prob);
   Result.Solver='fminsearch';

   % Call MATLAB OPTIM:s new simplex search-routine

   
   Result.SolverAlgorithm='fminsearch';

   if PriLev > 0
      fprintf('\n\nCall MATLAB Optimization TB routine fminsearch\n');
   end

   if exist('fminsearch','file')
      if any(isnan(Prob.x_0))
         Prob.x_0=[];
      end
      checkx0(Prob.x_0,Prob.x_L,Prob.x_U);
   
      if PriLev <= 0
         optstruct.Display= 'off';
      elseif PriLev == 0
         optstruct.Display= 'final';
      else
         optstruct.Display= 'iter';
      end
      optstruct.TolFun      = Prob.optParam.eps_f;
      optstruct.TolX        = Prob.optParam.eps_x;
      optstruct.MaxIter     = Prob.optParam.MaxIter;
      %optstruct.MaxFunEvals = 5*Prob.optParam.MaxIter;
      %optstruct.MaxFunEvals = Prob.optParam.MaxFunc;
      optstruct.MaxFunEvals =max(5*Prob.optParam.MaxIter,Prob.optParam.MaxFunc);

      if xnargin('fminsearch') == 6
         [x_k, Result.f_k, ExitFlag, out] = fminsearch(Prob.FUNCS.f, ...
               Prob.x_0, optstruct, Prob.x_L, Prob.x_U, Prob);
      else
         [x_k, Result.f_k, ExitFlag, out] = fminsearch('nlp_f', Prob.x_0, ...
               optstruct, Prob);
      end

      Result.Inform=ExitFlag;
      ExitFlag=max(-1,ExitFlag);
      switch ExitFlag
         case 0
           ExitFlag=1;  % Too many iterations
         case 1
           ExitFlag=0;  % OK
         case 2
           ExitFlag=0;  % OK
         case -1
           ExitFlag=3;  % Failure, but what type???
         case 10000
           ExitFlag=2;  % Unbounded ???
         case 10001
           ExitFlag=4;  % Infeasible ???
         case 10002
           ExitFlag=3;  % Rank problem ???
         case 10003
           ExitFlag=10; % Input errors ???
         otherwise
           ExitFlag=0;  % OK
      end

      Result.ExitFlag=ExitFlag;
      Result.x_k=x_k;
      Result.x_0 = Prob.x_0;
      if ~isempty(Prob.x_0)
         Result.f_0  = nlp_f(Prob.x_0, Prob);
      else
         Result.f_0 = NaN;
      end

      Result.FuncEv          = out.funcCount;
      Result.Iter            = out.iterations;
      Result.SolverAlgorithm = out.algorithm;

      % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
      %Result.xState=(x_k==Prob.x_L)+2*(x_k==Prob.x_U);
      Result=endSolve(Prob,Result);   
   else
      fprintf('Can not find quadprog.m.\n')
      fprintf('Have you got a license for Optimization Toolbox?\n')
      Result=[];
   end
elseif strcmpi(Solver,'LSQLIN') 

   % Call MATLAB OPTIM:s LS-routine

   Prob=ProbCheck(Prob,'LSQLIN',5);

   Prob=iniSolve(Prob,checkType('lls'),0,0);
   Result=ResultDef(Prob);
   Result.Solver='lsqlin';
   Result.SolverAlgorithm='LS algorithm in Optimization TB 2.x';

   if PriLev > 0
      fprintf('\n\nCall MATLAB Optimization TB routine lsqlin\n');
   end

   if any(isnan(Prob.x_0))
      Prob.x_0=[];
   end
   checkx0(Prob.x_0,Prob.x_L,Prob.x_U);
   if Prob.LineParam.LineAlg==1
      LineAlg='cubicpoly';
   else
      LineAlg='quadcubic';
   end
   if Prob.Solver.Alg==0
      % Levenberg-Marquardt
      alg='on';
      LargeScale='off';
   elseif Prob.Solver.Alg==1
      % Gauss Newton
      alg='off';
      LargeScale='off';
   else
      % Sparse
      alg='on';
      LargeScale='on';
   end
   if Prob.Solver.Method==0
      HessUpdate='bfgs';
   elseif Prob.Solver.Method==1
      HessUpdate='dfp';
   elseif Prob.Solver.Method==2
      HessUpdate='gillmurray';
   elseif Prob.Solver.Method==3
      HessUpdate='steepdesc';
   end
   if Prob.Solver.Alg==0
      optstruct.HessUpdate=HessUpdate;
      optstruct.LevenbergMarquardt=alg;
      optstruct.DiffMinChange=0.01*Prob.optParam.DiffInt;
      optstruct.DiffMaxChange=0.01;
   else
      optstruct.LargeScale= 'on';
   end

   if PriLev <= 0
      optstruct.Display= 'off';
   elseif PriLev == 0
      optstruct.Display= 'final';
   else
      optstruct.Display= 'iter';
   end
   optstruct.LargeScale     = LargeScale;
   optstruct.TolFun         = Prob.optParam.eps_f;
   optstruct.TolX           = Prob.optParam.eps_x;
   optstruct.MaxIter        = Prob.optParam.MaxIter;
   optstruct.TolCon         = Prob.optParam.cTol;
   optstruct.MaxFunEvals    = max(10000,Prob.optParam.MaxFunc);
   optstruct.Jacobian       = 'on';
   optstruct.LineSearchType = LineAlg;

   % Send info about functions to ls2_rJ
   Prob.FUNCSX.r   = Prob.FUNCS.r;
   Prob.FUNCSX.J   = Prob.FUNCS.J;
   Prob.FUNCSX.d2r = Prob.FUNCS.d2r;
   Prob.FUNCSX.rJ  = Prob.FUNCS.rJ;

   % Equalities
   Aeq = [];
   beq = [];
   if ~isempty(Prob.b_L)
      ix = find((Prob.b_L==Prob.b_U) & (abs(Prob.b_L) ~= Inf));
      if ~isempty(ix)
         Aeq = Prob.A(ix,:);
         beq = Prob.b_L(ix,:);
      end
   end
   if ~isempty(Prob.b_L)
      ix=find((Prob.b_L~=Prob.b_U) & (Prob.b_L ~= -Inf));
      if ~isempty(ix)
         A=-Prob.A(ix,:);
         b=-Prob.b_L(ix,:);
      else
         A=[];
         b=[];
      end
      ix=find((Prob.b_L~=Prob.b_U) & (Prob.b_U ~= Inf));
      if ~isempty(ix)
         A=[A; Prob.A(ix,:)];
         b=[b; Prob.b_U(ix,:)];
      end
   else
      ix = find(~isinf(Prob.b_U));
      A=Prob.A(ix,:);
      b=Prob.b_U(ix,:);
   end
      
   nA=xnargin('lsqlin');
   if nA > 10
      % Call TOMLAB version of lsqlin
      [x_k, Result.f_k, Result.r_k, ExitFlag, out, lam] = ...
            lsqlin(Prob.LS.C, Prob.LS.y,...
            A, b, Aeq, beq, Prob.x_L, Prob.x_U, Prob.x_0, ...
            optstruct, Prob);
   else
      [x_k, Result.f_k, Result.r_k, ExitFlag, out, lam] = ...
            lsqlin(Prob.LS.C, Prob.LS.y,...
            A, b, Aeq, beq, Prob.x_L, Prob.x_U, Prob.x_0, ...
            optstruct);
   end

   Result.Inform=ExitFlag;
   ExitFlag=max(-1,ExitFlag);
   switch ExitFlag
      case 0
        ExitFlag=1;  % Too many iterations
      case 1
        ExitFlag=0;  % OK
      case 2
        ExitFlag=0;  % OK
      case -1
        ExitFlag=3;  % Failure, but what type???
      case 10000
        ExitFlag=2;  % Unbounded ???
      case 10001
        ExitFlag=4;  % Infeasible ???
      case 10002
        ExitFlag=3;  % Rank problem ???
      case 10003
        ExitFlag=10; % Input errors ???
      otherwise
        ExitFlag=0;  % OK
   end

   Result.ExitFlag=ExitFlag;

   Result.f_k=0.5*Result.f_k; % OPT TB 2.0 returns r'*r, without 0.5 times.
   n = length(x_k);
   Result.v_k=[lam.lower;lam.upper];
   Result.x_k=x_k;
   Result.x_0 = Prob.x_0;
   Result.Iter=out.iterations;
   Result.FuncEv=out.iterations;
   Result.SolverAlgorithm=out.algorithm;

   % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
   Result.xState=(x_k(1:n)==Prob.x_L)+2*(x_k(1:n)==Prob.x_U);
   Result=endSolve(Prob,Result);
   
else
   Result = [];
   fprintf('Can not find the OPTIMIZATION TB 2.0 solver %s',Solver);
   fprintf('\n');
end
otxProb = [];

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

% MODIFICATION LOG:
%
% 990910  hkh  Add varargin to nonlinear solver calls.
% 011203  hkh  Add Result output when using Tomlab solvers 
% 011204  hkh  Handle different versions of fminsearch for Matlab 5.x and 6.x
% 030129  hkh  Remove varargin, always empty. Use global struct otxProb
% 030922  hkh  Set MaxFunEvals by use of Prob.optParam.MaxFunc
% 031129  hkh  n wrongly used after call (0 if x_0=[]), now n = Prob.N used
% 040415  hkh  Making correct calls to iniSolve
% 050422  hkh  Check if FUNCS.g,FUNCS.H is []; set options to fminunc/fmincon
% 050422  hkh  Check if FUNCS.J is []; set options to lsqnonlin
% 060814  med  FUNCS used for callbacks instead
% 060817  hkh  Added Prob.FUNCSX.rJ  = Prob.FUNCS.rJ; for LS solvers
% 080606  med  Cleaned
% 090804  med  Removed version check