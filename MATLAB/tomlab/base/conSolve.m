% conSolve.m:
%
% conSolve implements SQP algorithms for constrained minimization problems
%
% Main algorithm:
% A SQP algorithm by Schittkowski with Augmented Lagrangian merit function.
%
% conSolve also has an alternative Han-Powell SQP algorithm implemented
% (does not work as well as the Schittkowski algorithm)
%
% function Result = conSolve(Prob, varargin)
%
% Minimization problem:
%        min          f(x)
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%              c_L <= c(x) <= c_U
%
% Algorithms:
% conSolve runs different methods to obtain the gradient g and the Hessian H,
% dependent on the parameters Prob.Solver.Alg, Prob.NumDiff, Prob.ADObj
% and if a user supplied routine for the Hessian, stored in
% Prob.FUNCS.H is available.
%
% Solver.Alg=2 or 4 gives quasi-Newton BFGS as the Hessian approximation
%
% The table gives the different possibilities
%
% Solver.Alg == 0 gives the conSolve default algorithm choice.
%
% Schittkowski SQP:
% Solver.Alg  NumDiff   ADObj  isempty(FUNCS.H)  Hessian computation
%     0          0         0           0         Analytic Hessian
%     0        any       any         any         BFGS
%
%     1          0         0           0         Analytic Hessian
%     1          0         0           1         Numerical differences H
%     1         >0         0         any         Numerical differences g,H
%     1         <0         0         any         Numerical differences H
%     1        any        -1         any         Automatic differentiation
%     2          0         0         any         BFGS
%     2        ~=0         0         any         BFGS, numerical gradient g
%     2        any        -1         any         BFGS, automatic diff gradient
%
% Han-Powell SQP:
% Solver.Alg  NumDiff    ADObj isempty(FUNCS.H)  Hessian computation
%     3          0         0           0         Analytic Hessian
%     3          0         0           1         Numerical differences H
%     3         >0         0         any         Numerical differences g,H
%     3         <0         0         any         Numerical differences H
%     3        any         1         any         Automatic differentiation
%     4          0         0         any         BFGS
%     4        ~=0         0         any         BFGS, numerical gradient g
%     4        any         1         any         BFGS, automatic diff gradient
%
% INPUT PARAMETERS
%   Use conAssign.m (or probAssign) to initialize the Prob structure in the
%   TOMLAB Quick format.  Fields used in structure Prob:
%
%   x_0:     Starting point
%   x_L:     Lower bounds for x
%   x_U:     Upper bounds for x
%   b_L:     Lower bounds for linear constraints
%   b_U:     Upper bounds for linear constraints
%   A:       Linear constraint matrix
%   c_L:     Lower bounds for nonlinear constraints
%   c_U:     Upper bounds for nonlinear constraints
%   NumDiff  How to obtain derivatives (gradient, Hessian)
%   ConsDiff Differentiation method for the constraint Jacobian
%            0 = analytic, 1-5 different numerical methods
%   SolverQP Name of the solver used for QP subproblems. If empty,
%            picked from a list, best available with a license
%   f_Low    A lower bound on the optimal function value,
%            see LineParam.fLowBnd below. 
%            Used in convergence tests, f_k(x_k) <= f_Low. Only a feasible 
%            point x_k is accepted
%
%   FUNCS.f:   The routine to compute the function, as a string
%   FUNCS.g:   The routine to compute the gradient, as a string
%   FUNCS.H:   The routine to compute the Hessian, as a string
%   FUNCS.c:   The routine to evaluate the constraints, as a string
%   FUNCS.dc:  The routine to compute the gradient of the constraints
%
%   PriLevOpt Print level: 0 Silent, 1 Final result, 2 Each iteration
%   	      3 Line search info and Hessian
%
% optParam structure in Prob. Fields used:
%    bTol          Linear constraint violation convergence tolerance
%    cTol          Constraint violation convergence tolerance
%    eps_g         Gradient convergence tolerance
%    eps_x         Convergence tolerance in x
%    eps_Rank      Rank tolerance
%    IterPrint     Print short information each iteration
%    MaxIter       Maximal number of iterations
%    QN_InitMatrix Initial Quasi-Newton matrix, if not empty,
%                  Otherwise use identity matrix
%    size_x        Approximate size of optimal variable values
%    size_f        Approximate size of optimal f(x) value
%    xTol          Variable violation tolerance
%    wait          Pause after printout if true
%
%   LineParam Structure in Prob for line search parameters
%                See help LineSearch. One parameter is set in conAssign:
%   fLowBnd      A lower bound on the global optimum of f(x).
%                The user might also give lower bound estimate in Prob.f_Low
%                conSolve computes LineParam.fLowBnd as:
%                LineParam.fLowBnd = max(Prob.f_Low,Prob.LineParam.fLowBnd)
%
% Extra parameters:
%   VARARGIN : User defined parameters passed to f, c, g, dc and H
%
% OUTPUT PARAMETERS
% Result      Structure with results from optimization
%    x_k      Optimal point
%    v_k      Lagrange multipliers NOT USED
%    f_k      Function value at optimum
%    g_k      Gradient vector at optimum
%    H_k      Hessian value at optimum
%    x_0      Starting value vector
%    f_0      Function value at start
%    c_k      Constraint values at optimum
%    cJac     Constraint derivative values at optimum
%    xState   Variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
%    bState   Linear constraint: Inactive==0; On lower bound == 1;
%             On upper bound == 2; Equality == 3;
%    cState   Nonlinear constraint: Inactive==0; On lower bound == 1;
%             On upper bound == 2; Equality == 3;
%    Iter     Number of iterations
%    ExitFlag Flag giving exit status
%    ExitTest Text string giving ExitFlag and Inform information
%    Inform   Code telling type of convergence
%             1   Iteration points are close.
%             2   Small search direction.
%             3   Iteration points are close and small search direction.
%             4   Merit function gradient small.
%             5   Iteration points are close and merit function gradient
%                 small.
%             6   Small search direction and merit function gradient small.
%             7   Iteration points are close, small search direction and
%                 merit function gradient small
%             8   Small p and constraints satisfied.
%           101   Max no of iterations reached
%           102   Function value below given estimate.
%           103   Close iterations, but constraints not fulfilled.
%                 Too large penalty weights to be able to continue
%                 Problem is maybe infeasible?
%           104   Search direction 0 and infeasible constraints.
%                 The problem is very likely infeasible
%           105   Merit function is infinity.
%           106   Penalty weights too high
%
% Reference for the Schittkowski SQP:
%
% Klaus Schittkowski: On the convergence of a Sequential Quadratic Programming
% Method with an Augmented Lagrangian Line Search Function
% Systems Optimization Laboratory, Dept. of Operations Research, Stanford
% University, Stanford CA 94305, January 1982.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1995-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written June 5, 1995.   Last modified June 7, 2007.

function Result = conSolve(Prob, varargin)

%#function qp_f qp_g qp_H

if nargin < 1
   error('conSolve needs input structure Prob');
end

solvType=checkType('con');

Prob=ProbCheck(Prob,'conSolve',solvType);

global MAX_x % Max number of variables to print

[Alg, SolvAlg, B_k, merit, minpen]=NameAlg(Prob);
if Alg==2 | Alg==4
   Prob.optParam.ObjDerLev  = 1;
   Prob.optParam.ConsDerLev = 1;
   Prob = iniSolve(Prob,solvType,1,1);
else
   Prob.optParam.ObjDerLev  = 2;
   Prob.optParam.ConsDerLev = 1;
   Prob = iniSolve(Prob,solvType,2,1);
end

% If true (==1) plot line search problem at each iteration
plotLine=Prob.plotLine;

% Pick up input variables from Prob structure
optParam  = Prob.optParam;
LineParam = Prob.LineParam;
LineParam.fLowBnd = max(Prob.LineParam.fLowBnd,Prob.f_Low); 
fLow      = LineParam.fLowBnd;

[n, x_k, x_km1, xEqual, x_L, x_U, Prob] = BoundInit(Prob);

% Linear constraints

[mA, Ax, bEqual, b_L, b_U, A] = LinearConstr(Prob);

% Nonlinear constraints

[m, c_k, dc_k, cEqual, c_L, c_U] = NonlinConstr(Prob, varargin{:});

mTot=mA+m;

eps_g =  optParam.eps_g;    % Gradient convergence tolerance
epsRank= optParam.eps_Rank; % Rank test tolerance
eps_x =  optParam.eps_x;    % Convergence tolerance in x
MaxIter= optParam.MaxIter;  % Maximal number of iterations

size_x = optParam.size_x;   % Approximate size of optimal variable values 
size_f = optParam.size_f;   % Approximate size of optimal function value 
%size_c = optParam.size_c;   % Approximate size of optimal constraint value 

% If x in [x_L,x_L+xTol] or [x_U-xTol,x_U], fix x on bounds
xTol =  optParam.xTol;    
bTol =  optParam.bTol;    
cTol =  optParam.cTol;    

PriLev    = Prob.PriLevOpt;       % Print level
IterPrint = optParam.IterPrint;   % Print short information each iteration

if MaxIter <= 0, MaxIter=100*n; end
global xLast
if LineParam.LineAlg > 1
   CURV=1;
   LineParam.LineAlg=LineParam.LineAlg-2;
   xLast=x_km1;
else
   CURV=0;   % No curvi linear search
end

RESTART=1;  % RESTART == 1 now, no attempt to restart with bad penalties.

set = zeros(n,1);

Result=ResultDef(Prob);

rho_k=1;
rho_bar=100;                     %     rho_bar   > 1
delta_bar=0.9;                   % 0 < delta_bar < 1
%delta_k=0;                      % 0 <= delta_k <= 1  

Prob.Mode   = 2;
Prob.nState = 1;
f_k=nlp_f(x_k, Prob, varargin{:} );  % Function value
Prob.Mode   = 1;
g_k=nlp_g(x_k, Prob, varargin{:} );  % Gradient

if length(g_k)~=n
   fprintf('Number of variables %d. Length of gradient %d\n',...
            n,length(g_k));
   error('conSolve: Illegal length on gradient');
end

Result.x_0=x_k;
Result.f_0=f_k;

Bk_max=1E10*max(diag(B_k));

if Alg==1 | Alg==3
   H_k = nlp_H(x_k,Prob, varargin{:});% Hessian
   if isempty(H_k)
      % Switch to BFGS 
      Alg=Alg+2;
      SolvAlg=[SolvAlg, '!!! BFGS used'];
      B_k=speye(n,n);
   else
      B_k=H_k;
      if all(B_k(:) < eps), B_k=speye(n,n); end
   end
end
Prob.nState = 0;

if isempty(Prob.SolverQP)
   Convex = Alg == 2 | Alg == 4;
   if Convex
      if checkLicense('sqopt')
         SolverQP='sqopt';
      elseif checkLicense('lssol')
         SolverQP='lssol';
      elseif checkLicense('qpopt')
         SolverQP='qpopt';
      elseif checkLicense('bqpd')
         SolverQP='bqpd';
      else
         SolverQP='qld';
      end
   else
      if checkLicense('qpopt')
         SolverQP='qpopt';
      elseif checkLicense('bqpd')
         SolverQP='bqpd';
      else
         SolverQP='qpSolve';
      end
   end
   % SolverQP=GetSolver('qp',Prob.LargeScale,Convex);
   Prob.SolverQP = SolverQP;
else
   SolverQP=Prob.SolverQP;
end

ProbQP = [];

Result.Solver='conSolve';

Result.SolverAlgorithm=SolvAlg;

if Alg==3 | Alg==4  % Han-Powell
   %pen=0.001*ones(n+mTot,1);              % Penalty parameter
   pen=0.1*ones(n+mTot,1);              % Penalty parameter
else
   %pen=100*ones(n+mTot,1);              % Penalty parameter
   pen=ones(n+mTot,1);              % Penalty parameter
end

v_k=zeros(n+mTot,1);                % Estimate of Lagrange multipliers

bl=cat(1,x_L,b_L,c_L);
bu=cat(1,x_U,b_U,c_U);
zTol=[xTol*ones(n,1);bTol*ones(mA,1);cTol*ones(m,1)];
iE=cat(1,xEqual,bEqual,cEqual);
M=length(bl);
I=find(~iE);
E=find( iE);

% The field Prob.f_k etc. is used for communication with con_fm.m
Prob.f_k=f_k;
Prob.g_k=g_k;
Prob.c_k=c_k;
Prob.cJac=dc_k;
Prob.MERIT.Ax=Ax;
Prob.MERIT.pen=pen;
Prob.MERIT.merit=merit;
Prob.MERIT.bl=bl;
Prob.MERIT.bu=bu;
Prob.MERIT.E=E;
Prob.MERIT.I=I;

% Flag con_fm.m & con_gm.m that f_k, g_k, Ax, c_k & dc_k computed
Prob.MERIT.fceval=-1; 

Prob.Mode = 2;
if Alg==1 | Alg==2
   fM_k = con_fm([x_k;v_k],[],Prob,varargin{:});
   gM_k = con_gm([x_k;v_k],[],Prob,varargin{:});
else
   fM_k = con_fm(x_k,[],Prob,varargin{:});
   gM_k = con_gm(x_k,[],Prob,varargin{:});
end


maxpen=10000;
zBz=1;
p=Inf*ones(n,1);
global n_f
if IterPrint
   fprintf('Iteration Function         f(x)           |step|     ');
   fprintf('line search   |gradient|  Constraint\n');
   fprintf('           Count                                       ');
   fprintf('   step        (of meritf)  Violations\n');
end
p    =0;
gProj=0;
alpha=0;

Iter = 0;                       % Iteration index
while 1 % Main loop starts --------------------------------------------------
   zz=cat(1,x_k,Ax,c_k);
   h_k = norm(max(0,bl-zz),1) + norm(max(0,zz-bu),1);    
   if PriLev > 1
      fprintf('==================================================\n');
      fprintf('Iteration no: %4.0f  Function value %30.20f\n',Iter,f_k);
      if M > n
         fprintf('             Merit  Function value %30.20f\n',fM_k);
         fprintf('             f_k  + sum(|constr|)  %30.20f\n',f_k+h_k);
      end
   end
   if PriLev > 2
      xprint(x_k,'x_k:');
      xprinte(g_k,'g_k:');
      xprinte(gM_k,'gM_k');
      xprinte(c_k,'c_k:');
      xprinte(v_k,'v_k:');
   end
   if PriLev > 3
      xprinte(dc_k,'dc_k');
   end


   % SQ1. [Check termination criteria.]
% ---------------------------------------------------------------------------
%   CHECK CONVERGENCE HERE
   stop = 0;
   Inform = 0;
   flag = NaN;
   if (Alg==1 | Alg==2) & max(abs(x_k-x_km1)./max(abs(x_k),size_x)) <= eps_x   
      Inform=Inform+1;
      flag=0;
      stop = 1;
   end

   if zBz <= eps_x	% Convergence criteria 2
      flag = 0;
      Inform = Inform + 2;
      stop = 1;
   end
   %if max(abs(gM_k(idx)).* max(abs(x_k(idx)),size_x)) <= 
   if Alg==1 | Alg==2
      if max(abs(gM_k).* max(abs([x_k;v_k]),size_x)) <= ...
           eps_g * max(abs(f_k),size_f)
         flag=0;
         Inform = Inform + 4;
         stop = 1;
      end
   else
      if max(abs(gM_k).* max(abs(x_k),size_x)) <= ...
           eps_g * max(abs(f_k),size_f)
         flag=0;
         Inform = Inform + 4;
         stop = 1;
      end
   end

   if Iter >= MaxIter		% Stop criteria 1
      if Inform >0 & Inform < 100, flag=0; 
      else
          flag=1; Inform=101; 
      end
      stop = 1;
   end
   if f_k <= fLow      % Stop criteria 2
      if Inform >0 & Inform < 100, flag=0; 
      else 
          flag=1; Inform=102; 
      end
      stop = 1;
   end
   if isinf(fM_k)
      stop=1;  % We give up
      Inform=105;
      flag=1;
   end
   if (Inform < 100 | Inform==102) & stop
      % Check if nonlinear constraints are fulfilled
      
      cErr=max(0,max(c_k-c_U,c_L-c_k));
      ix=cErr > cTol*max(1,abs(c_k));
      ceq=sum(ix & cEqual);
      cineq=sum(ix & ~cEqual);
      ix=find(cErr > cTol*max(1,abs(c_k)));
      pen(n+mA+find(ix))=100*pen(n+mA+find(ix));

      if ceq > 0 & PriLev > 1
         fprintf('%d   equalities off more than cTol =%15.3e\n',ceq,cTol);
      end
      if cineq > 0 & PriLev > 1
         fprintf('%d inequalities off more than cTol =%15.3e\n',cineq,cTol);
      end
      % Check if linear constraints are fulfilled
      
      bErr=max(0,max(Ax-b_U,b_L-Ax));
      ix=bErr > bTol*max(1,abs(Ax));
      beq=sum(ix & bEqual);
      bineq=sum(ix & ~bEqual);
      pen(n+find(ix))=100*pen(n+find(ix));

      if beq > 0 & PriLev > 1
         fprintf('%d linear   equalities off > bTol =%15.3e\n',beq,bTol);
      end
      if bineq > 0 & PriLev > 1
         fprintf('%d linear inequalities off > bTol =%15.3e\n',bineq,bTol);
      end

      if (any([cineq,ceq,beq,bineq] > 0)) & ...
          (Iter < MaxIter & f_k >= fLow)
         % We have not minimum, but what to do???
         if any(pen > 1E18)
            stop=1;  % We give up
            Inform=103;
            flag=1;
         elseif all(p == 0)
            stop=1;  % We give up. Infeasible problem
            Inform=104;
            flag=1;
         elseif isinf(fM_k)
            stop=1;  % We give up
            Inform=105;
            flag=1;
         else
            stop=0; 
            if PriLev > 1
               fprintf('Equalities not fulfilled, Increase penalties\n');
            end
            minpen=2; % Emergency increase of penalties to speed convergence
            Prob.MERIT.pen=pen;

            Prob.f_k=f_k;
            Prob.g_k=g_k;
            Prob.c_k=c_k;
            Prob.cJac=dc_k;
            Prob.MERIT.Ax=Ax;
            Prob.MERIT.fceval=-1; % Flag to avoid recomputing
            Prob.Mode = 2;
            if Alg==1 | Alg==2
               fM_k = con_fm([x_k;v_k],[],Prob,varargin{:});
               gM_k = con_gm([x_k;v_k],[],Prob,varargin{:});
            else
               fM_k = con_fm(x_k,[],Prob,varargin{:});
               gM_k = con_gm(x_k,[],Prob,varargin{:});
            end
         end
      end

   end
   if ~stop & any(pen > 1E20) & ~RESTART
      RESTART=1;
      pen=10*ones(length(pen),1);
      rho_k=1;
      minpen=2;
   end
   if any(pen > 1E20) & RESTART
      stop=1;  % We give up
      Inform=106;
      flag=1;
   end 
   if IterPrint
      fprintf(' %5d  %7d   %20.17f %11.7e %10.6f %17.7e %10.5e\n', ...
              Iter, n_f, f_k, norm(p), alpha, norm(gM_k), h_k);
   end
   if stop            % Show final results
      if PriLev >= 1  
         switch Inform
            case 101
               disp('*** STOP 1! Max no of iterations reached ***');
            case 103
               disp('*** STOP 3! Too high penalty values ***');
            case 104
               disp('*** STOP 4! No progress. Infeasible problem? ***');
            case 105
               disp('*** STOP 4! Merit function is Inf ***');
            case 106
               disp('*** STOP 3! Too high penalty values ***');
            case 102
               disp('*** STOP 2! Function value below lower bound.');
               disp('*** Restart with lower Prob.f_Low ');
               disp('*** if minimum not reached ***');
            case 1
               disp('*** Convergence 1, Iteration points are close ***');
            case 2
               disp('*** Convergence 2, Small search direction ***');
            case 3
               disp('*** Convergence 1, Iteration points are close ***');
               disp('*** Convergence 2, Small search direction ***');
            case 4
               disp('*** Convergence 3, Merit function gradient small ***');
            case 5
               disp('*** Convergence 1, Iteration points are close ***');
               disp('*** Convergence 3, Merit function gradient small ***');
            case 6
               disp('*** Convergence 2, Small search direction ***');
               disp('*** Convergence 3, Merit function gradient small ***');
            case 7
               disp('*** Convergence 1, Iteration points are close ***');
               disp('*** Convergence 2, Small search direction ***');
               disp('*** Convergence 3, Merit function gradient small ***');
         end
      end
      if PriLev >= 5  
         fprintf('==================================================\n');
         fprintf('Iteration no: %4.0f  Function value %30.20f\n',Iter,f_k);
         xprint(x_k,'x_k:');
         xprinte(g_k,'g_k:');
         xprinte(gM_k,'gM_k');
         xprinte(c_k,'c_k:');
         xprinte(v_k,'v_k:');
         if Alg==1 | Alg==3
            PrintMatrix(B_k(1:min(n,MAX_x),1:min(n,MAX_x)),'The Hessian H_k:')
         end
      end
      Result.x_k=x_k;
      Result.v_k=v_k;

      Result.Iter=Iter;
      Result.ExitFlag=flag;
      Result.ExitText=ExitText(Inform);
      Result.Inform=Inform;
      
      Result.f_k=f_k;
      Result.g_k=g_k;
      Result.c_k=c_k;
      Result.cJac=dc_k;
      if Alg==1 | Alg==3
         Result.H_k=B_k;
      else
         Result.B_k=B_k;
      end
      % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
      Result.xState=double((x_k==x_L)+2*(x_k==x_U));
      if ~isempty(Prob.A)
         Result.bState=double((Ax<=b_L)+2*(Ax>=b_U));
      else
         Result.bState=[];
      end
      if ~isempty(c_k)
         Result.cState=double((c_k<=c_L)+2*(c_k>=c_U));
      else
         Result.cState=[];
      end

      Result=endSolve(Prob,Result);

      return
   end

% ---------------------------------------------------------------------------
%  END OF CONVERGENCE CHECK

   % Compute the Hessian if used
   if PriLev > 3 & (Alg==1 | Alg==3)
      PrintMatrix(B_k(1:min(n,MAX_x),1:min(n,MAX_x)),'The Hessian H_k:')
   end


   % SQ2. [Solve the quadratic programming subproblem.]

% ---------------------------------------------------------------------------
   Inform = 0;
   rho_0 = 1;
   [p, u_k, err, ProbQP, rho_max, Z]=QPsetup( B_k, x_k, x_L, x_U,...
 zTol, v_k, g_k, A, dc_k, c_k, bl, bu, Ax, rho_k, SolverQP, ProbQP, PriLev);

   fp_0=1;
   goon=1;
   while fp_0 > 0 & ~stop & goon
% ---------------------------------------------------------------------------
%     Some failure in the QP solution

      delta_k=0;
      pAccept=err==0;
      %while err > 0 & rho_k < rho_max
      while err == 4 & rho_k < rho_max
         ProbQP.QP.F(n+1,n+1)=rho_k;

         ResultQP = tomRunMini(SolverQP,ProbQP);

         p=ResultQP.x_k;
         u_k=ResultQP.v_k;
         u_k=[u_k(1:n);u_k(n+2:M+1)];
         err=ResultQP.ExitFlag;

         delta_k=p(n+1);
         if delta_bar <= delta_k & (err==0)
            % Check if rank problems
            [Q, R, EE, ppRank] = ComputeQR(Z(:,1:n), epsRank);
            if ppRank < min(size(Z,1),n) & any(p(1:n)~=0)
               % Accept trial p. 
               pAccept=1;
               delta_k=0.5; % Must set below 1, to avoid division with 0
            else
               pAccept=0;
            end
         elseif delta_bar > delta_k & (err==0)
            pAccept=1;
         else
            pAccept=0;
         end
         if pAccept
            p=p(1:n);

         else
            rho_k=rho_bar*rho_k;
         end
      end
%err
%rho_k
%pAccept
%pause
      %if rho_k > rho_max
      %   rho_k=Bk_max;
      %end
      
      if ~pAccept & (Alg==1 | Alg==2) % Take gradient step in Merit function
         if PriLev > 1
            fprintf('conSolve: Take gradient step in Merit function.\n');
            fprintf('QP-ErrorFlag=%d. Iter =%d, delta_k=%15.5e.',...
                    err,Iter,delta_k);
            fprintf('\n');
         end
         % Reinitialize penalty weights after gradient step if BIG
         if any(pen > 1E15)
            pen = 100*ones(length(pen),1);
         end
         % p=-B_k\gM_k(1:n);

         [Q, R, EE, pRank] = ComputeQR(B_k, epsRank);
         p = tomsol(6, Q, R, EE, pRank, -gM_k(1:n));

         if PriLev > 1
            fprintf('QR with Pivot: Dim %d Rank %3.0f\n',size(B_k,2),pRank);
         end

         % Check if variables on bounds, and avoid the wrong direction
         zL=abs(x_k-x_L) < xTol*max(1,abs(x_k));
         zU=abs(x_k-x_U) < xTol*max(1,abs(x_k));
         z0=zL|zU;
         p(find(z0))=0;
         %u_k=v_k-gM_k(n+1:length(gM_k));
         u_k=v_k+sign(v_k)*norm(v_k)/max(1E-5,norm(p));
         % Add search steps in Lagrange Params for the vars on bounds???
         %u_k(find(zL))=1;
         %u_k(find(zU))=-1;
         p_v=u_k-v_k;
         delta_k=0;
         % Restrict the length of the step to a reasonable level
         xNorm=max(1,norm(x_k));
         pNorm=norm(p);
         if pNorm > 10*xNorm, p=xNorm*p/pNorm; end
      elseif ~pAccept & (Alg==3 | Alg==4)
         z=cat(1,Ax,c_k);
         inAct=~(bl(n+1:M)==bu(n+1:M)) & ~(v_k(n+1:M) > 0 | ...
           min(z-bl(n+1:M),bu(n+1:M)-z) <= zTol(n+1:M).*max(1,abs(z)));
         xAct=min(x_k-x_L,x_U-x_k) <= xTol*max(1,abs(x_k));
         Active=iE | [xAct;~inAct];
         dz=[speye(n,n);A;dc_k]';
         %p=-dz(:,find(Active))\g_k;
         % These lines are not necessary
         % i=find((bu-[x_k;z]) < ([x_k;z]-bl));
         % Active(i)=-Active(i);
         % dz(:,i)=-dz(:,i); % Change sign when upper bounds active

         % Compute projected gradient, should get descent
         ixx=find(Active);
         [Qg, Rg, Eg, pRank] = ComputeQR(dz(:,ixx), epsRank);

         %NOTp2=-(speye(n,n)-dz(:,find(Active))*pinv(dz(:,find(Active))))*g_k;
         % Compute projected gradient
         p = -(g_k - Qg(:,1:pRank)*Qg(:,1:pRank)'*g_k);
         clear Qg Rg Eg

         % Restrict the length of the step to a reasonable level
         xNorm=max(1,norm(x_k));
         pNorm=norm(p);
         if pNorm > 10*xNorm, p=xNorm*p/pNorm; end
         fp_0=gM_k'*p;
         if fp_0 >= 0
            fp_0=-fp_0;  % Just change sign to get descent
            p=-p;
         end
      else

         if isempty(A)
            z=[x_k+p;c_k];
         else
            z=[x_k+p;A*(x_k+p);c_k];
         end
         epsL=1E-14;
         test1=all(abs(p) <= eps_x*max(1,abs(x_k)));
         %test2=sum(abs(c_k(1:me))) + sum(abs(min(0,c_k(me+1:M)))) <= cTol;
         test2=sum(abs(z(E)-bl(E))) + sum(abs(min(0,min(z(I)-bl(I),...
                  bu(I)-z(I))))) <= min(zTol);
         %Xtest1=norm(p) <= eps_x;
         %X%test2=sum(abs(c_k(1:me))) + sum(abs(min(0,c_k(me+1:M)))) <= cTol;
         %Xtest2=sum(abs(z(E)-bl(E))) + sum(abs(min(0,min(z(I)-bl(I),...
         %X         bu(I)-z(I))))) <= zTol(I).*max(1,abs(z(I)));
         % HKH??? Must some sign dependence on constraint be used here?
         test3=0&all(v_k(I)+epsL >= 0);
         %if norm(p) <= eps_x & ...
         %   (sum(abs(c_k(1:me))) + sum(abs(min(0,c_k(me+1:M)))) <= cTol) & ...
         %   all(v_k(me+1:M) >= -E-14) 
         if test1 & test2 & test3
            if PriLev >= 1
               fprintf('*** Convergence 4, Small p & constraints satisfied\n');
            end
            p_v=zeros(M,1);
            flag=0;
            Inform = Inform + 8;
            stop = 1;
         else
            p_v=u_k-v_k;

            if ~test1
               pBp=p'*B_k*p;
               if pBp==0 | delta_k==1
                  z=2*length(pen);
               else
                  z=2*length(pen)/((1-delta_k)*pBp);
               end
            else
               z=2*length(pen);
            end
            if Alg==1 | Alg==2
               sig_k=1-(1-1./sqrt(pen)).^(Iter+2);
               %sig_k=min(1,(Iter+2)./sqrt(pen)); %Old simplified Schittkowski

               for i=1:length(pen) % Update penalty parameters
                   % HKH special emergency update of penalties

                   % Restrict u-v-square, could blow away
                   uv= min(100*sig_k(i)*pen(i), z*(u_k(i)-v_k(i))^2);
                   pen(i)=min(maxpen*pen(i),...
                              max([minpen*pen(i),sig_k(i)*pen(i),uv]));

                   % Works bad if square blows
                   %pen(i)=min(maxpen*pen(i),...
                   %           max([minpen*pen(i),sig_k(i)*pen(i), ...
                   %                z*(u_k(i)-v_k(i))^2]));

                   % Original Schittkowski
                   %pen(i)=max(sig_k(i)*pen(i),z*(u_k(i)-v_k(i))^2);
               end
            else
               v_kNorm=norm(v_k,inf);
               for i=1:length(pen) % Update penalty parameters
                   pen(i)=max([minpen*pen(i),v_kNorm]);
                   % pen(i)=max([minpen*pen(i),abs(v_k(i)),1]);
               end
            end

            Prob.MERIT.pen=pen;

            Prob.f_k=f_k;
            Prob.g_k=g_k;
            Prob.c_k=c_k;
            Prob.cJac=dc_k;
            Prob.MERIT.Ax=Ax;
            Prob.MERIT.fceval=-1; % Flag to avoid recomputing
            Prob.Mode = 2;
            if Alg==1 | Alg==2
               fM_k = con_fm([x_k;v_k],[],Prob,varargin{:});
               gM_k = con_gm([x_k;v_k],[],Prob,varargin{:});
            else
               fM_k = con_fm(x_k,[],Prob,varargin{:});
               gM_k = con_gm(x_k,[],Prob,varargin{:});
            end
         end
      end
      if ~stop
         if Alg==1 | Alg==2
            fp_0=gM_k'*[p(:);p_v(:)];
         else
            fp_0=gM_k'*p;
         end
         if fp_0 > 100*eps & ~pAccept
            % Failure in gradient step
            if PriLev > 1
               disp('Failure in gradient step')
            end
            goon=0;
         elseif abs(fp_0) <= 100*eps & err==0
            % Do nothing, accept this
            goon=0;
            if norm(p) < eps_x
               stop=0;
            end
             
         else
            if fp_0 >= 0, rho_k=rho_bar*rho_k; end
            err=10; 
         end
      end
   end
% ---------------------------------------------------------------------------

   if stop			
      Result.x_k=x_k;
      Result.v_k=v_k;

      Result.Iter=Iter;
      Result.ExitFlag=flag;
      Result.ExitText=ExitText(Inform);
      Result.Inform=Inform;
      
      Result.f_k=f_k;
      Result.g_k=g_k;
      Result.c_k=c_k;
      Result.cJac=dc_k;
      % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
      Result.xState=double((x_k==x_L)+2*(x_k==x_U));
      if ~isempty(Prob.A)
         Result.bState=double((Ax<=b_L)+2*(Ax>=b_U));
      else
         Result.bState=[];
      end
      if ~isempty(c_k)
         Result.cState=double((c_k<=c_L)+2*(c_k>=c_U));
      else
         Result.cState=[];
      end

      Result=endSolve(Prob,Result);

      return % Terminate
   elseif delta_k < 1 & ~all(p==0)
      % Adaptive change of initial value of rho_k, used in case of infeasible QP
      z=p'*[speye(n,n),A',dc_k']*u_k;
      % CHANGED rho_k=max(rho_k,10*z^2/((1-delta_k)^2*p'*B_k*p)); 
      pBp=p'*B_k*p;
      if pBp==0 | delta_k==1
         rho_k=rho_0;
      else
         rho_k=max(rho_0,10*z^2/((1-delta_k)^2*pBp)); 
      end
   end
% ---------------------------------------------------------------------------
 

   % SQ3. [Update the estimate of the solution]
   % Update parameters
   f_km1 = f_k;
   g_km1 = g_k;
   c_km1 = c_k;
   dc_km1 = dc_k;
   fM_km1 = fM_k;
   gM_km1 = gM_k;

   alpha_1=LineParam.InitStepLength;
   if PriLev > 2
      fprintf('p:');
      xprinte(p);
      xprint(p,'p:  ');
      if PriLev > 2
         if alpha_1 < 0.99
            fprintf('Line search start %10.6f\n',alpha_1)
         end
      end
   end
   alphaMax = 1E20;
   p(xEqual)=0;
   for i=1:n
       if p(i) > xTol*max(1,abs(x_k))
          alphaMax = min(alphaMax,(x_U(i)-x_k(i))/p(i));
       else
          if p(i) < -xTol*max(1,abs(x_k))
             alphaMax = min(alphaMax,(x_L(i)-x_k(i))/p(i));
          end
       end
   end
   if PriLev > 2
      fprintf('Alpha max  %10.6f\n',alphaMax)
   end
   Prob.f_k=f_k;
   Prob.g_k=g_k;
   Prob.c_k=c_k;
   Prob.cJac=dc_k;

   alphaMax=min(3,alphaMax);

   Prob.MERIT.fceval=1; % Flag that function values NOT sent

   LineParam.InitStepLength = alpha_1;
   if Alg==1 | Alg==2
      LineResult = LineSearch('con_fm','con_gm',...
          [x_k;v_k], [p;p_v], fM_km1, gM_km1, LineParam, alphaMax, ...
           3, PriLev-3, Prob, varargin{:});

      alpha=LineResult.alpha;
      if plotLine
         Prob.Mode = 0;
         LinePlot('con_fm', [x_k;v_k], [p;p_v], fM_km1, gM_km1,...
                   LineParam, alphaMax, 3, alpha, Prob, varargin{:});
      end
   else
      LineResult = LineSearch('con_fm','con_gm',...
          x_k, p, fM_km1, gM_km1, LineParam, alphaMax, ...
           3, PriLev-3, Prob, varargin{:});

      alpha=LineResult.alpha;
      if plotLine
         Prob.Mode = 0;
         LinePlot('con_fm', x_k, p, fM_km1, gM_km1,...
                   LineParam, alphaMax, 3, alpha, Prob, varargin{:});
      end
   end

   fM_k    = LineResult.f_alpha;
   gM_k    = LineResult.g_alpha;
   alphaVec= LineResult.alphaVec;

   if alpha > 1E-14 
      % Pick up values computed from con_fm and con_gm
      f_k = LineResult.f_k;  % Function value

      if ~isempty(LineResult.g_k)
         g_k = LineResult.g_k;  % Gradient
      end
      if ~isempty(LineResult.Ax)
         Ax = LineResult.Ax(:);  
      end
      if ~isempty(LineResult.c_k) & ~isempty(LineResult.cJac) 
         c_k = LineResult.c_k;  % Constraints 
         dc_k= LineResult.cJac; % Constraint Jacobian
      end
   end


   if PriLev > 2
      fprintf('Line search in search direction p,')
      fprintf(' Step length alpha = %12.6f\n',alpha);
   end
% ---------------------------------------------------------------------------

   if optParam.wait, pause, end % Please wait!
 
   % Update parameters
   x_km1 = x_k;
   x_k = x_k + alpha*p;
   if ~isempty(A) 
      Ax=A*x_k;
      Prob.MERIT.Ax=Ax;
   end
   if Alg==1 | Alg==2
      v_k = v_k + alpha*p_v;
   else
      v_k = u_k;
   end

   if Alg==1 | Alg==3
      B_k = nlp_H(x_k,Prob, varargin{:}); % Hessian
      if all(B_k(:) < eps), B_k=speye(n,n); end
   else
      [B_k, idx]=SafeUpdate(alpha,p,eps_x,eps_g,B_k,g_k,g_km1,PriLev,Iter);
      B_k = 0.5*(B_k+B_k');
   end

   Iter = Iter+1;
   if max(pen) > 1E10
      minpen=1000;
   end
end % Main loop ends --------------------------------------------------------


% =========================================
function [B_k, idx]=SafeUpdate(alpha,p,eps_x,eps_g,B_k,g_k,g_km1,PriLev,Iter)
% =========================================
eps1000 = 1000*eps;    % 1000 * machine precision

global MAX_x % Max number of variables to print

   % Safeguarded BFGS update of approximate Hessian
   z=alpha*p; % Check if we have converged (or short alpha step)
   n=length(p);
   idx=[1:n];  % Maybe restrict which variables to update?
   if norm(z) > eps_x
      y=g_k(idx)-g_km1(idx);
      if norm(y) > eps_g
         Bz=B_k(idx,idx)*z;
         zBz=z'*Bz;
         zy=z'*y;
         if zy < 0.2 * zBz
            theta=0.8*zBz / (zBz - zy);
            w=theta * y + (1-theta)*Bz;
            zw=z'*w;
         else 
            w=y;
            zw=zy;
         end
         if PriLev > 2
            fprintf('BFGS Hessian update iter %3.0f\n',Iter);
            fprintf('zw %10.6e zBz %10.6e \n',zw,zBz)
         end
         if zw < 1E-13 | zBz < 1E-13
            if PriLev > 1
               fprintf('BFGS Hessian update step dangerous!\n');
               fprintf('max(w)/zw %20.10e\n',max(w)/zw);
               fprintf('min(w)/zw %20.10e\n',min(w)/zw);
               fprintf('max(Bz)/zBz %20.10e\n',max(Bz)/zBz);
               fprintf('min(Bz)/zBz %20.10e\n',min(Bz)/zBz);
            end
         end
         if zw == 0  
            if zBz ~=0
               B_k(idx,idx)=B_k(idx,idx)-(Bz/zBz)*Bz';
            end
         elseif  zBz == 0
            B_k(idx,idx)=B_k(idx,idx)+(w/zw)*w';
         else
            B_k(idx,idx)=B_k(idx,idx)+(w/zw)*w'-(Bz/zBz)*Bz';
         end
         if PriLev > 3
            fprintf('B_k matrix after BFGS update\n');
            PrintMatrix(B_k(1:min(n,MAX_x),1:min(n,MAX_x)),'B_k:')
         end
         if (PriLev > 2 & n < 50) | PriLev > 10
            fprintf('B_k matrix eigenvalues:\n');
            if issparse(B_k)
               xprint(full(eigs(B_k)),'eig:')
            else
               xprint(eig(B_k),'eig:')
            end
         end
      end
   end
% =========================================
function [Alg, SolvAlg, B_k, merit, minpen]=NameAlg(Prob)
% =========================================

% Algorithm. 0= Default alogrithm, Schittkowski SQP, analytic, or BFGS
%            1= Schittkowski SQP, with analytic Hessian
%               Schittkowski SQP, with numerical Hessian, Hessian not supplied
%            2= Schittkowski SQP with BFGS update,
%            3= Han-Powell SQP, with numerical or analytical Hessian
%            4= Han-Powell SQP with BFGS update.

Alg=max(0,Prob.Solver.Alg);

if isempty(Alg), Alg=0; end
if Alg > 4, Alg=0; end

if Alg==1 & Prob.NumDiff==0 & (~isempty(Prob.FUNCS.H) | any(Prob.ADObj == -1))
   SolvAlg='Schittkowski SQP algorithm, analytic Hessian';
   merit=0; % merit == 0 ==> Augmented Lagrangian
   minpen=0;   % HKH special. Emergency increase of the penalities
elseif Alg==1 
   SolvAlg='Schittkowski SQP algorithm, Numerical Hessian';
   merit=0; % merit == 0 ==> Augmented Lagrangian
   minpen=0;   % HKH special. Emergency increase of the penalities
elseif Alg==2
   SolvAlg='Schittkowski SQP algorithm with BFGS update';
   merit=0; % merit == 0 ==> Augmented Lagrangian
   minpen=0;   % HKH special. Emergency increase of the penalities
elseif Alg==3 & Prob.NumDiff==0 & ...
       (~isempty(Prob.FUNCS.H) | any(Prob.ADObj == -1))
   SolvAlg='Han-Powell SQP algorithm';
   merit=1; % merit == 1 ==> Non-differentiable L1 merit function
   minpen=2;
elseif Alg==3 
   SolvAlg='Han-Powell SQP algorithm with Numerical Hessian';
   merit=1; % merit == 1 ==> Non-differentiable L1 merit function
   minpen=2;
elseif Alg==4 
   SolvAlg='Han-Powell SQP algorithm with BFGS update';
   merit=1; % merit == 1 ==> Non-differentiable L1 merit function
   minpen=2;
else
   % Default algorithm
   if Prob.NumDiff==0 & (~isempty(Prob.FUNCS.H) | any(Prob.ADObj == -1))
      SolvAlg='Schittkowski SQP algorithm, analytic Hessian';
      merit=0; % merit == 0 ==> Augmented Lagrangian
      minpen=0;   % HKH special. Emergency increase of the penalities
      Alg = 1;
   %elseif isempty(Prob.FUNCS.g) | Prob.NumDiff > 0
   else
      SolvAlg='Schittkowski SQP algorithm with BFGS update';
      merit=0; % merit == 0 ==> Augmented Lagrangian
      minpen=0;   % HKH special. Emergency increase of the penalities
      Alg = 2;
   end
end

if Alg==2 | Alg==4
   % Find initial Quasi-Newton matrix
   n = length(Prob.x_0);
   if ~isempty(Prob.optParam.QN_InitMatrix)
      B_k = Prob.optParam.QN_InitMatrix;   % B = User given matrix
      if size(B_k,1)~=n | size(B_k,2)~=n 
         disp('Illegal dimensions on initial Quasi Newton matrix');
         disp('conSolve is using identity matrix instead');
   
         B_k = speye(n,n);                 % B = Use identity matrix
      end
   else
      B_k = speye(n,n); % B_k = Start Quasi-Newton with identity matrix
   end
else
   B_k=[];
end

% =========================================
function [p, u_k, err, ProbQP, rho_max, Z]=QPsetup(B_k, x_k, x_L, x_U,...
   zTol, v_k, g_k, A, dc_k, c_k, bl, bu, Ax, rho_k, SolverQP, ProbQP, PriLev)
% =========================================

   % SQ2. [Solve the quadratic programming subproblem.]

% ---------------------------------------------------------------------------


stop    = 0;
rho_0   = 1;
Z       = [];
rho_max = [];

n       = length(x_k);
m       = length(c_k);
mA      = size(A,1);
M       = length(bl);

if isempty(A)
   AA = dc_k; 
elseif isempty(dc_k)
   AA = A; 
else
   AA = [A;dc_k]; 
end
if M > n
   z   = cat(1,Ax,c_k);
   b_L = bl(n+1:M)-z;
   b_U = bu(n+1:M)-z;
else
   b_L = [];
   b_U = [];
end

if isempty(ProbQP)
   % Setup Structure used in QP calls
   % ProbQP = CreateProbQP(Prob, 'qp', max(500,3*n), PriLev-3, optParam.wait);
   % Setup structure used in QP call
   ProbQP = qpAssign(B_k,g_k,AA,b_L,b_U,x_L-x_k,x_U-x_k,zeros(n,1),'QP');
   ProbQP.optParam         = optParamDef(SolverQP,2,n,0, M-n);
   ProbQP.optParam.MaxIter = max([500,3*n,ProbQP.optParam.MaxIter]);
else
   ProbQP.QP.F = B_k;
   ProbQP.QP.c = g_k;
   ProbQP.N    = n;
   ProbQP.x_L  = x_L-x_k;
   ProbQP.x_U  = x_U-x_k;
   ProbQP.x_0  = zeros(n,1);
   ProbQP.b_L  = b_L;
   ProbQP.b_U  = b_U;
   ProbQP.A    = AA;
   ProbQP.mLin = M-n;
end

ResultQP = tomRunMini(SolverQP,ProbQP);
p   = ResultQP.x_k;
err = ResultQP.ExitFlag;
if norm(p)==0   % Total failure, try B_k == I, gradient step
   ProbQP.QP.F = speye(n,n);
   ProbQP.Name = 'GradQP';
   ResultQP    = tomRunMini(SolverQP,ProbQP);
   p           = ResultQP.x_k;
   err         = ResultQP.ExitFlag;
   ProbQP.Name = 'QP';
end
u_k = ResultQP.v_k;
if err==1
   % Max iterations reached. Maybe the solution is useful
   fp_0=g_k'*p;
   if fp_0 < -zTol(1)
      % Accept this solution
      err=0;
   else
   end
end


%%% THINK ABOUT THIS!!! No linear constraints
%%ProbQP.b_L=[];
%%ProbQP.b_U=[];
%%% No constraints
%%err=0;
%%p=-pinv(full(B_k))*g_k;
%%step=Inf;
%%for i=1:n
%%    if p(i) < 0
%%       step=min(step,(ProbQP.x_L(i)-x_k(i))/p(i));
%%    elseif p(i) > 0
%%       step=min(step,(ProbQP.x_U(i)-x_k(i))/p(i));
%%    end
%%end
%%step=min(step,1);
%%u_k=zeros(n,1);
%%p=step*p;
%%gg=B_k*p+g_k;
%%for i=1:n
%%    if p(i) <= ProbQP.x_L(i)
%%       u_k(i)=gg(i);
%%    elseif p(i) >= ProbQP.x_U(i)
%%       u_k(i)=-gg(i);
%%    end
%%end



if err == 4
   if PriLev > 1
      fprintf('FAILURE TO SOLVE QP, err =%3.0f. Expand QP\n',err)
   end
   % HKH NO POINT IN COMPUTING QP solution once again
   %ResultQP = tomSolve(SolverQP,ProbQP,0);
   Bk_max   = 1E10*max(diag(B_k));
   rho_max  = min(1E15,max(200,1E10*Bk_max));
   % WHY COMPUTE ??? k=n+size(Ax,1);
   inAct=~(bl(n+1:M)==bu(n+1:M)) &...
     ~(v_k(n+1:M) > 0 | min(z-bl(n+1:M),bu(n+1:M)-z) <= ...
     zTol(n+1:M).*max(1,abs(z)));

   z(find(inAct))=0;

   % Expand problem
   ProbQP.Name = 'Expanded QP';
   ProbQP.QP.F = [ProbQP.QP.F, zeros(n,1);zeros(1,n),0];
   ProbQP.QP.c = [g_k;0];
   ProbQP.N    = n+1;
   Z           = [ProbQP.A,zeros(m+mA,1)];

   % Lower bound -Inf  is still -Inf. Correct upper bound is already set
   iLow             = find(isinf(bl(n+1:M)) & ~inAct);
   %ProbQP.b_U(iLow)=bu(n+iLow)-z;
   % Subtract [b_U,c_U] for new last column in matrix
   z(iLow)          = z(iLow)-bu(n+iLow);
   ProbQP.b_L(iLow) = -z(iLow);

   % Upper bound  Inf  is still  Inf. Correct lower bound is already set
   iUpp             = find(isinf(bu(n+1:M)) & ~inAct);
   %ProbQP.b_L(iUpp)=bl(n+iUpp)-z;
   % Subtract [b_L,c_L] for new last column in matrix
   z(iUpp)          = z(iUpp)-bl(n+iUpp);
   ProbQP.b_L(iUpp) = -z(iUpp);

   % Equalities. Correct bounds. Must correct for new last column
   iEQ              = find(bl(n+1:M)==bu(n+1:M) & ~inAct);
   z(iEQ)           = z(iEQ)-bl(n+iEQ);
   ProbQP.b_L(iEQ)  = -z(iEQ);
   ProbQP.b_U(iEQ)  = -z(iEQ);

   % Active inequalities, where no bound isInf
   % Lower bound active
   inE1=find((z-bl(n+1:M) <  bu(n+1:M)-z) & ~inAct & ~isinf(bl(n+1:M))...
             & ~isinf(bu(n+1:M)) & ~(bl(n+1:M)==bu(n+1:M)));
   if ~isempty(inE1)
      z(inE1)=z(inE1)-bl(n+inE1);
      ProbQP.b_L(inE1) = -z(inE1);
      ProbQP.b_U(inE1) = Inf;
   end
 
   inE2=find((z-bl(n+1:M) >= bu(n+1:M)-z) & ~inAct & ~isinf(bl(n+1:M))...
             & ~isinf(bu(n+1:M)) & ~(bl(n+1:M)==bu(n+1:M)));
   % Upper bound active
   if ~isempty(inE2)
      z(inE2)          = z(inE2)-bu(n+inE2);
      ProbQP.b_U(inE2) = -z(inE2);
      ProbQP.b_L(inE2) = -Inf;
   end

   % Add the extra column for the extra delta variable
   % Add extra constraint rows when no bound isInf 
   ProbQP.A            = [[ProbQP.A,-z];Z(inE1,:);Z(inE2,:)];
   ProbQP.mLin         = size(ProbQP.A,1);

   %if ~isempty(inE1)
   %   if ~isempty(inE2)
   %      ProbQP.A=[[ProbQP.A,z];[Z(inE1,:),0];[Z(inE2,:),0]];
   %   else
   %      ProbQP.A=[[ProbQP.A,z];[Z(inE1,:),0]];
   %   end
   %elseif ~isempty(inE2)
   %   ProbQP.A=[[ProbQP.A,z];[Z(inE2,:),0]];
   %else
   %   ProbQP.A=[ProbQP.A,z];
   %end

   % Add bounds for the NOT active constraints when no bound isInf 
   ProbQP.b_L=[ProbQP.b_L;...
              -Inf*ones(length(inE1),1);...
               bl(n+inE2);
               ];
   ProbQP.b_U=[ProbQP.b_U;...
               bu(n+inE1);...
               Inf*ones(length(inE2),1);...
               ];
   % Add extra column for the extra delta variable

   ProbQP.x_L = [ProbQP.x_L;0];
   ProbQP.x_U = [ProbQP.x_U;1];
   ProbQP.x_0 = [ProbQP.x_0;1];

   %zU=z+bu(n+1:M);
   %iLow=find(~isinf(bl(n+1:M)));
   % Must make special treatment of upper bounds for inequalities
   % For active equations with lower bound -Inf, use upper bound equations
   %iUpp=find((~(bl(n+1:M)==bu(n+1:M)) & ~inAct &...
   %          ~isinf(bu(n+1:M))) | isinf(bl(n+1:M));
   %ProbQP.A=[ProbQP.A([iLow,iUpp],:),[zL(iLow);zU(iUpp)]];
   % Now also using linear equations. Maybe only use nonlinear ???????
   %if ~isempty(iUpp)
   %   % Add extra equations for active inequalities with valid upper bounds
   %   ProbQP.b_U=[ProbQP.b_U;ProbQP.b_U(iUpp,:)];
   %   % Eliminate upper bounds for corresponding equations
   %   ProbQP.b_U(iUpp)=Inf;
   %   % Add extra rows in A for active inequalities 
   %   ProbQP.A=[ProbQP.A;[ProbQP.A(iUpp,1:M-n),(z(iUpp)+bu(n+iUpp))]];
   %   % Dummy lower bounds for  active inequalities 
   %   ProbQP.b_L=[ProbQP.b_L;-Inf*ones(length(iUpp),1)];
   %end
   %ProbQP.A=[A;dc_k]; 
   %ProbQP.b_L=bl(n+1:M)-[Ax;c_k];
   %ProbQP.b_U=bu(n+1:M)-[Ax;c_k];
else
   rho_max=100*rho_k;
end

% ------------------------------
function Text = ExitText(Inform)
% ------------------------------

switch  Inform
   case 1
     Text = 'Iteration points are close';
   case 2
     Text = 'Small search direction';
   case 3
     Text = 'Iterations close, small search direction';
   case 4
     Text = 'Merit function gradient small';
   case 5
     Text = 'Iterations close, merit function gradient small';
   case 6
     Text = 'Small search direction, merit function gradient small';
   case 7
     Text = str2mat('Iterations close, small search direction' ...
                   ,'merit function gradient small');
   case 8
     Text = 'Small search direction, constraints satisfied';
   case 32
     Text = 'Local minimum with all variables on bounds';
   case 101
     Text = 'Maximal number of iterations reached';
   case 102
     Text = str2mat('Function value below lower bound fLow' ...
                   ,'Restart with lower fLow if minimum not reached');
   case 103
     Text = str2mat('Close iterations, but constraints not fulfilled' ...
                   ,'Too large penalty weights to be able to continue' ...
                   ,'Problem is maybe infeasible');
   case 104
     Text = str2mat('Zero search direction and infeasible constraints' ...
                   ,'Problem is very likely infeasible');
   case 105
     Text = 'Merit function is Infinity';
   case 106
     Text = 'Too high penalty values';
end
if Inform < 100
   Text=str2mat('Optimal solution found',Text);
end

% MODIFICATION LOG:
%
% 981001  hkh  Changed to structure Result and Prob
% 981005  mbk  Changes in comments concerning type of convergence.
% 981005  mbk  ucDef changed to conDef on line 118.
% 981006  hkh  New call to LineSearch. Using LinePlot. Cleanup.
% 981006  hkh  Better convergence tests.
% 981013  hkh  Added call to iniSolve and endSolve
% 981016  hkh  LP problems made B_k singular. Must check if y small before
%              update
% 981017  hkh  Deleted the gradient as argument to LinePlot
%              Changed the logic in input/output to merit function gradient
% 981019  hkh  Added Han-Powell algorithm to conSolve
% 981020  hkh  Add Projected gradient step as last resort in Han-Powell
% 981021  hkh  Check if qpopt is present, otherwise use qpSolve.
% 981026  hkh  Restrict gradient step to reasonable level. Improve printing.
% 981026  hkh  Solver is name of solver. SolverAlgorithm is description
%              Use optParam.Method to make choice of QPOPT or qpSolve
%              Change Alg, Alg=0 Schittkowski, Alg=1, Han-Powell
% 981027  hkh  Use Prob.Solver.Alg to define Alg
% 981028  hkh  Put f_0 = f(x_0) as field in Result
% 981110  hkh  Change LineResult.dc_k to cJac
% 981110  hkh  Change to MAX_x,MAX_c,MAX_r
% 981112  hkh  Add check on empty x_0, and try finding some reasonable x_0
% 981117  hkh  Use QR to compute p if PANIC
% 981128  hkh  Change some printing levels for output
% 981203  hkh  Must check min(n,MAX) in output
% 990706  hkh  Expand algorithms to handle numerical Hessian and BFGS.
% 990824  hkh  Use general routine SolverQP for QP subproblems
% 000916  hkh  Adding text for ExitText
% 001014  hkh  Use Prob.ConsDiff for diff. method for the constraint Jacobian
% 001015  hkh  Reset penalty if emergency gradient step, improves convergence
% 001101  hkh  Use qld as default QP solver for v2.1 and v3.0 /MINI
% 020409  hkh  Use sparse eye (speye) everywhere. Important if big zero H.
% 020409  hkh  Set Prob.Mode and Prob.nState before function/gradient calls
% 020701  hkh  Must set B_k full in p=-pinv(full(B_k))*g_k;
% 020820  ago  Change for Matlab 6.5 logical handling
% 030211  hkh  Change comment for eps_absf, to clearify
% 031128  med  Modifications for MAD, new ADObj and ADCons
% 031129  hkh  Call NameAlg before iniSolve. Compute ObjDerLev for selected alg
% 031201  hkh  Revise algorithm selection for use of MAD and ADMAT
% 040111  hkh  Change call to inisolve
% 040125  hkh  Set fields mLin in ProbQP before call to tomSolve
% 040414  hkh  Correct call to iniSolve regarding levels of derivatives
% 040728  med  Pragmas added for MATLAB Compiler
% 041124  hkh  Safe guard pen update for (u_k(i)-v_k(i))^2 blowing up
% 041125  hkh  Formulate QP to solve unconstrained step, previous strategy
%              not foolproof, still left in code
% 041126  hkh  Safe guard rho_max <= 1E15, avoiding crashes in qpopt
% 041126  hkh  Only expand QP problem if ExitFlag == 4, not ExitFlag > 0
% 051216  med  Help updated, empty print removed
% 060814  med  FUNCS used for callbacks instead
% 060818  hkh  Set LineParam.fLowBnd = max(Prob.LineParam.fLowBnd,Prob.f_Low); 
% 060818  hkh  Use Prob.f_Low instead of optParam.eps_absf for target f test
% 061212  med  ADMAT removed
% 070221  med  Help updated
% 070712  frhe Make B_k symmetric, or CPLEX complains.
% 070907  hkh  SolverQP picked from list, 1st with license, avoid GetSolver 
% 080607  hkh  Use tomSolve, not tomSolve, use qpAssign, not CreateProbQP
% 080607  med  Switched to tomRunMini
