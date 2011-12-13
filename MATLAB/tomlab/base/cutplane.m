% cutplane.m:
%
% function Result = cutplane(Prob)
%
% Cutting plane algorithm for Mixed-Integer Programming (MIP)
% using LP relaxation (Formulated as min-IP)
%
% Solving MIP on the LP standard form, with additional upper bounds on x:
%
%        min    c' * x.  x in R^n, A x = b_U, x >= 0, x <= x_U
%         x
%
% Any of the x could be set as integer valued
%
% Note that the TOMLAB routine cpTransf.m can be used to rewrite
% an LP problem on another form to this standard form.
%
%
% INPUT PARAMETERS
% Fields in Prob:
%   c:      The vector c in c'x.
%   A:      The linear constraint matrix
%   b_L:    Lower bounds on linear constraints. NOT USED. Assumed to be == b_U
%   b_U:    The upper bounds for the linear constraints
%   x_L:    Lower bounds on x. NOT USED. Assumed to be 0.
%   x_U:    Upper bounds on x. If empty assumed to be Inf.
%
%   x_0:    Starting point x (If EMPTY, lpSimplex solves Phase I LP to find a
%           feasible point.
% ---------------------------------------
% QP.B:     Active set B_0 at start.
% ---------------------------------------
%           1  = Include variable x(i) is in basic set.
%           0  = Variable x(i) is set on its lower bound
%           -1 = Variable x(i) is set on its upper bound. NOT USED
%           If EMPTY, lpSimplex finds active set.
%
% ---------------------------------------
% MIP       Structure in Prob, Prob.MIP.
% ---------------------------------------
%           Defines integer optimization parameters. Fields used:
%   IntVars:  
%           If empty, all variables are assumed non-integer 
%           If islogical(IntVars) (=all elements are 0/1), then
%           1 = integer variable, 0 = continuous variable.
%           If any element >1, IntVars is the indices for integer variables
%
% Prob.Solver.Alg     Not used
% Prob.Solver.Method  Rule to select new variables:
%           = 0 Minimum Reduced Cost, sort variables increasing. (DEFAULT)
%           = 1 Bland's Anti-cycling Rule
%           = 2 Minimum Reduced Cost, Dantzig's rule
% PriLev    Print level in cutplane
% SolverLP  Name of the solver used for initial LP subproblem. If empty,
%           picked from a list, best available with a license
% SolverDLP Name of the solver used for dual LP subproblems. If empty,
%           SolverLP is used.
%
% ----------------------------------------------------------------------
%
% Fields used in Prob.optParam, in lpSimplex and DualSolve
% PriLev    Printing level (in lpSimplex and DualSolve):
%           =0 No output; >0 Convergence results;
%           >1 Output every iteration  >2 Output each step in the simplex alg
% MaxIter   Maximal number of iterations. max(10*dim(x),100) is DEFAULT.
% wait      Wait flag, pause each iteration if set true
% eps_f     Tolerance used to test convergence on the reduced costs
% eps_Rank  Rank tolerance
% xTol      Tolerance to judge if x-values are close
%
% OUTPUT PARAMETERS
% Structure Result. Fields used:
%   Iter     Number of iterations
%   ExitFlag Exit flag
%            == 0  => OK
%            == 1  => Maximal number of iterations reached. No bfs found.
%            == 4  => No feasible point found with Phase1 Simplex
%            Otherwise  the same flags as the dual LP solver returns
%   ExitTest Text string giving ExitFlag and Inform information
%   Inform   If ExitFlag > 0, Inform=ExitFlag.
%   x_k      Solution
%   v_k      Lagrange parameters. Constraints + lower + upper bounds
%   QP.B     B  Optimal set. B(i)==1, include variable x(i) in basic set.
%            sum(B==1)==length(b)  holds. See QP.B as input.
%   QP.y     Dual parameters y (also part of v_k)
%   p_dx     Search steps in x
%   alphaV   Step lengths for each search step
%   f_k      Function value c'*x
%   g_k      Gradient c
%   Solver   cutplane
%   SolverAlgorithm  Description of method used
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1994-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.2.0$
% Written Nov 28, 1994.   Last modified Jun 7, 2008.

function ResultMIP = cutplane(Prob)

if nargin < 1
   error('cutplane needs input structure Prob');
end

solvType=checkType('mip');

Prob=ProbCheck(Prob,'cutplane',solvType);

Prob = iniSolveMini(Prob);

A   = sparse(Prob.A);
b   = Prob.b_U(:);
c   = Prob.QP.c(:);

[m, n] = size(A);

x = Prob.x_0(:);  
if any(isinf(x) | isnan(x)), x=[]; end

B = Prob.QP.B(:);

x_L = zeros(n,1);
x_U = Prob.x_U(:);
if isempty(x_U), x_U=Inf*ones(n,1); end
% Safe-guard starting point
x   = max(x_L(1:n),min(x,x_U(1:n)));
b_L = b;

Prob.N   = n;
Prob.x_L = x_L;
Prob.x_U = x_U;
Prob.b_L = b;

optParam=Prob.optParam;

wait     = optParam.wait;     % Pause after printout if true
epsRank  = optParam.eps_Rank; % Rank test tolerance
MaxIter  = optParam.MaxIter;  % Maximal number of iterations
eps_f    = optParam.eps_f;    % Reduced cost convergence tolerance
xTol     = optParam.xTol;
IterPrint= optParam.IterPrint;

if isempty(xTol), xTol = 100*eps; end

PriLev    = Prob.PriLev;      % Print level
if isempty(PriLev), PriLev=2; end
PriLevOpt = Prob.PriLevOpt;  % Print level in lpSimplex and DualSolve
%PriLev=2
%PriLevOpt=2

% Integer variables
IntVars  = DefPar(Prob.MIP,'IntVars',[]);

% Logical vector for integers
IV = false(n,1);

if isempty(IntVars)
   % No binary variables B or integer variables of type I
elseif any(IntVars==0) | all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > n)
      error('cutplane: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars = find(IV);

nIV=length(IntVars);

eps_1 = 1E-8;  % Level when a variable is considered as an integer value

x_IP_max=[];
B_IP_max=[];

if ~isfield(Prob.MIP,'VarWeight')
   Prob.MIP.VarWeight = [];
   VarWeight=[];
else
   VarWeight=Prob.MIP.VarWeight(:);
end

ResultMIP                 = ResultDef(Prob);
ResultMIP.Solver          = 'cutplane';
ResultMIP.Prob            = Prob;
ResultMIP.SolverAlgorithm = 'Cutting plane with Gomory MIP-cuts.';

% Node selection method
Alg=Prob.Solver.Alg;

if isempty(Alg), Alg=0; end

Prob.Solver.Alg=Alg;

%if Alg==2
%   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ' Breadth First.'];
%elseif Alg==1
%   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ' Depth First.'];
%else
%   Alg=0;
%   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ...
%         ' Depth First, then Breadth.'];
%end

% Index selection rule in DualSolve and lpSimplex

Method=Prob.Solver.Method;

if isempty(Method), Method=0; end

Prob.Solver.Method=Method;

if Method==1
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ... 
         ' (Bland''s rule).'];
elseif Method==2
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ...
         ' (Minimum reduced cost. Dantzig''s rule).'];
else
   Method=0;
   ResultMIP.SolverAlgorithm=[ResultMIP.SolverAlgorithm ...
         ' (Minimum reduced cost).'];
end

if PriLev > 0
   fprintf('Run %s',ResultMIP.SolverAlgorithm);
   fprintf('\n');
end

IndexRule = Prob.Solver.Method;

ResultMIP.f_0 = [];
ResultMIP.x_k  = [];
ResultMIP.f_k  = [];
ResultMIP.v_k  = [];
ResultMIP.QP.B = [];
ResultMIP.QP.y = [];

if isempty(Prob.SolverLP)
   if checkLicense('minos')
      SolverLP='lp-minos';
   elseif checkLicense('bqpd')
      SolverLP='bqpd';
   elseif checkLicense('xa')
      SolverLP='xa';
   else
      SolverLP='milpSolve';
   end
   % SolverLP=GetSolver('lp',Prob.LargeScale);
   Prob.SolverLP = SolverLP;
else
   SolverLP=Prob.SolverLP;
end

if isempty(Prob.SolverDLP)
   SolverDLP = SolverLP;
   %if checkLicense('minos')
   %   SolverDLP='minos';
   %elseif checkLicense('bqpd')
   %   SolverDLP='bqpd';
   %elseif checkLicense('xa')
   %   SolverDLP='xa';
   %else
   %   SolverDLP='milpSolve';
   %end
   % SolverDLP=GetSolver('dlp',Prob.LargeScale);
   Prob.SolverDLP = SolverDLP;
else
   SolverDLP=Prob.SolverDLP;
end

% Setup Structure used in LP calls
ProbLP  = CreateProbQP(Prob, 'lp', max(2000,3*n), PriLev-3, optParam.wait);

ProbLP.b_U  = b;
ProbLP.b_L  = b; 
ProbLP.A    = A;
ProbLP.mLin = size(A,1);
ProbLP.QP.c = c;
ProbLP.x_L  = x_L; 
ProbLP.x_U  = x_U; 
ProbLP.P    = 0; 

%ProbLP.QP.B = B;
%ProbLP.x_0   =[];
%ProbLP.y    = [];
%ProbLP.QP.DualLimit = fIPMax; % Try stop dual iterations early

% Setup Structure used in dual LP calls
ProbDLP = CreateProbQP(Prob, 'dlp', max(1000,3*n), PriLev-3, optParam.wait);
ProbDLP.PriLevOpt   = PriLevOpt;

% Initial step
if PriLev > 2 & PriLevOpt >= 1
   disp('=+=+=+=+=+=+=+=')
   fprintf('=== cutplane:    LP relaxation. Call lpSimplex.\n')
   disp('=+=+=+=+=+=+=+=')
end

% Demand min.cost
if PriLev ==2 & PriLevOpt >=2
   ProbLP.PriLevOpt   = PriLevOpt-1;
   PriLevOpt=PriLevOpt-1;
else
   ProbLP.PriLevOpt   = PriLevOpt;
end

ResultLP = tomRunMini(SolverLP,ProbLP);

ExitFlag = ResultLP.ExitFlag;
x        = ResultLP.x_k;
v_k      = ResultLP.v_k;
if isempty(v_k)
   y     = [];
else
   y     = v_k(n+1:n+m);
end
B        = ResultLP.QP.B;
f_k      = ResultLP.f_k;
Iter     = ResultLP.Iter;


if ExitFlag > 0
   if PriLev >= 2
      disp('No Phase II solution found to LP relaxation')
   end
   ResultMIP.x_k  = x(1:n);
   ResultMIP.f_k  = f_k;
   ResultMIP.v_k  = v_k;
   ResultMIP.QP.B = B;
   ResultMIP.QP.y = y;
   ResultMIP.ExitFlag=4;
   ResultMIP.ExitText=ExitText(4);
   ResultMIP = endSolveMini(Prob,ResultMIP);
   return;
end

iB = find(B==0);
iN = find(B==1);
iU = find(B==-1);


cTx=c'*x;
cTx0=cTx;

if PriLev > 2
   disp('Phase II solution x:')
   xprint(x,'x:')
end
if PriLev > 1
   fprintf('cutplane: cTx =%40.20f\n',cTx)
end


n0=n;
m0=m;

ixU = find(~isinf(x_U));
mU  = length(ixU);

% Remove slack bounds

Apos=all( A' >= 0)';

% Tighten bounds

if any(Apos)
   Ukeep=ones(mU,1);
   for i=1:n
       ix=find(Apos & A(:,i) > 0);
       if ~isempty(ix)
          z=floor(b(ix)./A(ix,i)+eps_1);
          U=min(z);
          if U <= x_U(i)
             % Bound is implicit
             Ukeep(i)=0;
             %disp('Bound is implicit')
          end
      end
   end
   if any(Ukeep==0)
      ixU = ixU(Ukeep==1);
      mU  = length(ixU);
   end
end

if mU > 0
   % mU slack variables for the upper bound linear constraints
   L = zeros(1,n);
   L(IntVars)=1;
   uL=L(ixU);
   iuL = find(uL==1);
   % Skip forcing slack as integer variables
   % IntVars = [IntVars(:);n+iuL(:)];
   % nIV=length(IntVars);

   x_U(ixU(iuL))=floor(x_U(ixU(iuL)) + eps_1*max(1,abs(x_U(ixU(iuL)))));

   x_L = [x_L;zeros(mU,1)];
   x_U = [x_U;x_U(ixU)];

   if ~isempty(VarWeight)
      maxW = max(VarWeight(~isinf(VarWeight)));
      VarWeight = [VarWeight;maxW*ones(mU,1)];
   end
     
   % mU extra linear constraints for upper bounds
   x = [x;x_U(ixU)-x(ixU)];
   c = [c;zeros(mU,1)];
   B = [B(1:n0);B(ixU)==-1];
   B(B(1:n0)==-1)=0;  % Upper bound variables put into basis
   b   = [b;x_U(ixU)];

   A = sparse([A zeros(m,mU); ...
       full(sparse(1:mU,ixU,ones(mU,1),mU,n)) eye(mU,mU)]);

   n = n + mU;
   m = m + mU;

   %B = [B(1:n0);zeros(mU,1)];
   %B(n0+find(abs(x_U(n0+1:n0+mU)-x(n0+1:n0+mU)) < ...
   %         xTol.*max(1,abs(x(n0+1:n0+mU))))) = -1;
   %B(n0+find(abs(x_L(n0+1:n0+mU)-x(n0+1:n0+mU)) < ...
   %         xTol.*max(1,abs(x(n0+1:n0+mU))))) = 1;
end
x_U0=Inf*ones(n,1);
ix = x_L==x_U;
x_U0(ix)=x_U(ix);

ProbDLP.y    = [];
%ProbDLP.x_U  = x_U; 
ProbDLP.x_U  = Inf*ones(n,1); 

Iidx = IntVars(B(IntVars)==0);
x_I  = floor(x(Iidx) + eps_1*max(1,abs(x(Iidx))));
x_r  = max(0,x(Iidx) - x_I);
nReal = full(sum(x_r > eps_1*max(1,abs(x(Iidx)))));

L       = zeros(n,1);
L(IntVars) = 1;
nn=n;

Improve = 0;
OldC = 0;
Iter = 0;

while Iter <= MaxIter
   Iter = Iter + 1;
   ix   = find(x_L==x_U0);
   if ~isempty(ix)
      B(ix)=2;
      x(ix)=x_L(ix);
   end
   iN   = find(B==1);
   iB   = find(B<=0);

   mB   = length(iB);
   if mB < m
      % disp('MUST INCREASE RANK. FATAL ERROR');
      % Compute the QR-decomposition

      [Q, R, E, pRank] = ComputeQR(A(:,iB), epsRank);
      % Generate a full rank basis, if possible
      pRank0=pRank;
      [B, Q, R, E, pRank, DelVar] = FullRank(B, A, Q, R, E, pRank, m, epsRank);
      if PriLev > 1 & pRank0 < pRank
         fprintf('FullRank: Increased rank from %d to %d.',pRank0,pRank);
         fprintf(' Rows m = %d.\n',m);
      end
      if ~isempty(DelVar)
         x(DelVar)=0;
      end

      iN   = find(B==1);
      iB   = find(B<=0);
   end

%if 0
%xprinti(B,'B:')
%xprinti(iB,'iB:');
%xprinte(x(iB),'xiB');
%pause
%end
   if PriLev > 1
      fprintf('\n--- cutplane iter =%7.0d',Iter)
      fprintf('. f_k %17.10f. ',f_k)
      fprintf('Vars %d. Cons %d. ',length(B),mB)
      fprintf('Noninteger %d of %d\n',nReal,nIV)
   end

   % Check if solution is integer
   idx = find(B(IntVars)==0);
   Iidx  = IntVars(idx);


   if ~isempty(Iidx)
      x_I   = floor(x(Iidx) + eps_1*max(1,abs(x(Iidx))));
      x_r   = max(0,x(Iidx)-x_I);
      nReal = full(sum(x_r > eps_1*max(1,abs(x(Iidx)))));
      if isempty(VarWeight)
         % Variables with frac.part closest to 0.5
         %[best_frac xBest] = min(abs(x_r-0.5));
         r=abs(x_r-0.5);
         [rBest iBest]=sort(r);
         irBest=find(rBest < 0.5-eps_1);
      else
         r=VarWeight(Iidx);
         r(x_r < eps_1*max(1,x_r))=Inf;
         [rBest iBest]=sort(r);
         irBest=find(~isinf(rBest));
      end
      if ~isempty(irBest)
         NewC = 1;
         iCand = iBest(irBest);
         Cand  = Iidx(iCand);
%fprintf('l-iCand %d, l-OldC %d',length(irBest),length(OldC))
         for i = 1:length(iCand)
             j = find(Cand(i) == OldC);
             if isempty(j)
                NewC = i;
                OldC = [OldC,Cand(i)];
                break;
             end
         end
         if NewC == 1 
            OldC = Cand(1); 
         end

         %if length(Cand) > 1 & OldC == Cand(1)
         %   % Switch candidates
         %   t = iCand(1);
         %   iCand(1) = iCand(2);
         %   iCand(2) = t;
         %   t = Cand(1);
         %   Cand(1) = Cand(2);
         %   Cand(2) = t;
         %end
         

         if PriLev > 1 & ~isempty(irBest)
            fprintf('Fraction close to .5 = %13.10f',x_r(iCand(1)))
            fprintf(' Base# %4.0f. Var# %4.0f\n',iCand(1),Cand(1))
            fprintf('Picked = %13.10f',x_r(iCand(NewC)))
            fprintf(' Base# %4.0f. Var# %4.0f\n',iCand(NewC),Cand(NewC))
            if PriLev >2, xprint(x(Iidx),'x:'); end
         end
         if PriLev > 2
            xprinti(Cand,'Var:');
            xprinte(x_r(iCand),'x_r');
         end
      end
   else
      irBest=[];
   end
   if IterPrint & PriLev < 2
      fprintf('It %4.0d',Iter)
      fprintf('. f %17.16f ',f_k)
      fprintf('Vrs %d Cns %d ',length(B),mB)
      fprintf('Nonint %d/%d ',nReal,nIV)
      fprintf('Red %15.7e ',Improve)
      if ~isempty(irBest) 
         fprintf('ClFrac %6.4f ',x_r(iCand(1)))
         fprintf('B# %2.0f V# %2.0f ',iCand(1),Cand(1))
         fprintf('Pick %6.4f ',x_r(iCand(NewC)))
         fprintf('B# %2.0f V# %2.0f ',iCand(NewC),Cand(NewC))
      end
      fprintf('\n')
   end
   if isempty(irBest) 
      ResultMIP.x_k  = x(1:n0);
      ResultMIP.f_k  = cTx;
      ResultMIP.v_k  = v_k;
      ResultMIP.QP.B = B_IP_max;
      ResultMIP.QP.y = y;
      ResultMIP.ExitFlag=0;
      ResultMIP.ExitText=ExitText(0);

      ResultMIP.QP.A   = A;
      ResultMIP.QP.b_U = b;
      ResultMIP.QP.b_L = b;
      ResultMIP.QP.x_L = x_L;
      ResultMIP.QP.x_U = x_U;
      ResultMIP.QP.c   = c;

      if PriLev > 0
         fprintf('\n--- Cutting Plane converged!!! ITER = %5.0f\n',Iter)
         fprintf('\n    Optimal Objective function =  %30.16f\n',cTx);
         %if 0  % HKH DEBUG output
         %   bTy=b'*y;
         %   fprintf('\n\nMin LP c''*x =     %30.16f\n',cTx);
         %   fprintf('Max LP b''*y =     %30.16f\n',bTy);
         %   % cTx always equal to b'*y
         %   fprintf('Duality gap =     %30.16f\n',bTy-cTx);
         %   c_hat = zeros(n, 1);
         %   iN = find(B==1);
         %   c_hat(iN) = c(iN) - A(:, iN)' * y;
         %   if PriLev > 1
         %      fprintf('Reduced costs:\n');
         %      xprinte(c_hat(iN),'c_hat:');
         %   end
         %end
      end
      if PriLev > 1
         xprint(B,'B:',' %4.0f');
         xprinti(x(1:n0),'x:')
      end
      ResultMIP = endSolveMini(Prob,ResultMIP);
      return;
   end
 
   % Gomory mixed-integer cut for variable Cand(NewC) 
   % (Base variable iCand(NewC))
   ixC = iCand(NewC);
   xC  = Cand(NewC);
   %if NewC > 1 & length(OldC) > 5
   %   OldC = OldC(Cand(1) ~= OldC);
   %end
   f0  =  x_r(ixC);
   % ----- Solve the linear equation  A_B * C_N = A_N

   % Find variable in basis B
   qHat=find(xC==iB);
   
   e_p = zeros(length(iB),1);
   e_p(qHat) = 1;
   % ----- Solve the linear equation  A_B' * u = e_qHat. Row qHat in A_B^(-1)

   % QR - factorization already computed: [Q R E] = qr(full(A(:,iB)));

   %Q     = ResultLP.QP.Q;
   %R     = ResultLP.QP.R;
   %E     = ResultLP.QP.E;
   %pRank = ResultLP.QP.pRank;
   %u = SolveQRT(Q, R, E, pRank, e_p);

   % Compute the QR-decomposition

   [Q, R, E, pRank] = ComputeQR(A(:,iB)', epsRank);

   if pRank < length(iB)
      % disp('MUST INCREASE RANK AFTER QR. FATAL ERROR');
      % Compute the QR-decomposition

      [Q, R, E, pRank] = ComputeQR(A(:,iB), epsRank);
      % Generate a full rank basis, if possible
      pRank0=pRank;
      [B, Q, R, E, pRank, DelVar] = FullRank(B, A, Q, R, E, pRank, ...
             length(iB), epsRank);
      if PriLev > 1 & pRank0 < pRank
      %if PriLev > -1 
         fprintf('FullRank: Increased rank from %d to %d.',pRank0,pRank);
         fprintf(' Rows m = %d.\n',m);
      end
      if ~isempty(DelVar)
         x(DelVar)=0;
      end
      if pRank < length(iB)
         [P,jE]=find(E);
         ix =iB(P(pRank+1:length(iB)));
         B(ix)=1;
         x(ix)=0;
         m=pRank;
      end

      iN   = find(B==1);
      iB   = find(B<=0);
      [Q, R, E, pRank] = ComputeQR(A(:,iB)', epsRank);

      % Find variable in basis B
      qHat=find(xC==iB);
   
      e_p = zeros(length(iB),1);
      e_p(qHat) = 1;
   end


   u = tomsol(6, Q, R, E, pRank, e_p);

   %u2 = A(:,iB)' \ e_p

   c_N = u' * A(:, iN);

%AI=inv(A(:,iB));
%CN=AI(qHat,:)*A(:,iN);
%xprinte(CN,'CN');
%xprinte(c_N,'cn');
%sum(abs(CN-c_N))

%   if Iter==1
%      iB1=find(B==0);
%      AI=inv(A(:,iB1));
%   else
%      %iU = find(B==-1)-n0;
%      %mU =length(iU);
%      %if mU > 0
%      %   AA = [A  zeros(m0,mU); eye(mU,n0),eye(mU,mU)];
%      %end
%      AI=inv(A(:,iB));
%   end
%%u3=AI(qHat,:)
%pause
  
%   c_N0 = AI(qHat, :) * A(:, iN);
%   sum(abs(AI(qHat,:)'-u))

   a_i = zeros(1,n);
Williams=0;

if Williams
   b =  [b;1-f0];   % Create new constraint, rhs best fraction
   a_i(iN) = max(0,-c_N - floor(-c_N + eps_1));

   % Not integer + Not in base B. Reset to reduced cost. 
   ix = find(L==0 & B==1);

   a_i(ix) = -c_N(L(B==1)==0);

   % Not integer + Not in base B + Reduced cost < 0. Multiply by f0q.
   ix = find(L==0 & B==1 & a_i(:) < -eps_1);
   a_i(ix) = -(1-f0)/f0 * a_i(ix);
   a_i(a_i < eps_1)=0;
else
   b = [b; f0];        % Create new constraint, rhs best fraction
   f0q = f0 /(1-f0);
   % Compute fj for both integer and non-integer variables
   a_i(iN) = max(0,c_N - floor(c_N + eps_1));
   % These coefficients are only OK if fj <= f0 and c_N(Real vars) > 0
   % Compute the others:

   % integer + Not in base B + fraction > rhs fraction f0
   ix = find(L==1 & B==1 & a_i(:) > f0+eps_1);

   a_i(ix) = f0q*(1-a_i(ix));

   % Not integer + Not in base B. Reset to reduced cost. 
   ix = find(L==0 & B==1);

   a_i(ix) = c_N(L(B==1)==0);

   % Not integer + Not in base B + Reduced cost < 0. Multiply by f0q.
   ix = find(L==0 & B==1 & a_i(:) < -eps_1);
   a_i(ix) = -f0q * a_i(ix);

   a_i(a_i < eps_1)=0;
if 1
   %ix = find(a_i(1:n0-m0) > 0)';
   ix = find(a_i(1:n0) > 0)';
   if ~isempty(ix)
      z = floor(f0./a_i(1,ix)'+eps_1*max(1,x(ix)));
      iz = find(z < x_U(ix));
      if ~isempty(iz)
         % disp('Tighter bounds')
         % [ix(iz) x(ix(iz)) z(iz) x_U(ix(iz))]
         x_U(ix(iz))=z(iz);
         x_U0(ix(iz))=z(iz);
         ix=find(x_U0==x_L);
         if ~isempty(ix)
            x_U0(ix)=x_L(ix);
         end
         %pause
      end
   else
      ix = find(a_i(1:nn) > 0)';
      % xprinte(a_i,'a_i');
   end
end

end

if sum(a_i < 0) > 0
   disp('error. Some terms are negative')
   sum(a_i < 0)
   ix=find(a_i <0)
   a_i(ix)
   pause
end

if length(iB)~=pRank
   disp('ERROR - problem with rank')
   pRank
   length(iB)
   length(iN)
   n
   pause
end
   % Add new column
   A = [A,zeros(m,1); a_i,-1];
   B = [B;0];
   L = [L;0];
   c = [c;0];
   n = n+1;
   m = m+1;
   %x_U(n)   = b(m);
   %x_U0(n)  = b(m);
   x_U(n)   = Inf;
   x_U0(n)  = Inf;
   x_L(n)   = 0;
   x        = full([x;b(m)]);
   x(iN)    = 0;
   x(x<0)   = 0;
   ix       = x > x_U;
   if ~isempty(ix)
      x(ix) = x_U(ix);
   end

   if PriLev > 2
      xprint(b,'b:');
   end
   if PriLev > 1 & PriLevOpt >= 1
      disp('=+=+=+=+=+=+=+=')
      disp('=== cutplane:     Call dual LP solver:')
      disp('=+=+=+=+=+=+=+=')
   end

   %if Iter > 6
   %   ProbDLP.PriLevOpt   = 4;
   %end
   %ProbDLP.PriLevOpt   = 2;
   %ProbDLP.optParam.wait     = 1;
   %ProbDLP.PriLevOpt   = 3;

   ProbDLP.b_U   = b;
   ProbDLP.b_L   = b; % DualSolve must know all eq:s are equalities
   ProbDLP.A     = A;
   ProbDLP.mLin  = size(A,1);
   ProbDLP.QP.c  = c;
   ProbDLP.QP.B  = B;
   ProbDLP.x_0   = x;
   ProbDLP.x_L   = x_L; 
   ProbDLP.x_U   = x_U0;
   %ProbDLP.x_U  = [Inf*ones(nn,1); x_U(nn+1:n)];
   %ProbDLP.x_U  = x_U;
   %ProbDLP.QP.DualLimit = fIP_max; % Try stop dual iterations early
   ProbDLP.P  = Iter; 
   ProbDLP.N  = length(x); 

   %if mod(Iter,50)==0
   %   xxU       = ProbDLP.x_U;
   %   ProbDLP.x_U = x_U;
   %   Probx     = preSolve(ProbDLP);
   %   fprintf('Sum nan %d ',sum(isnan(Probx.b_L)));
   %   fprintf('Sum fixed %d',sum(Probx.x_L==Probx.x_U));
   %   ProbDLP.x_U = xxU;
   %end

   ResultLP = tomRunMini(SolverDLP,ProbDLP);

   ExitFlag=ResultLP.ExitFlag;
   
   x   = ResultLP.x_k;
   %y   = v_k(n+1:length(v_k));
   fDP = ResultLP.f_k;
   B   = ResultLP.QP.B;

   %B(abs(x-x_U) < xTol*max(1,abs(x)))=-1;
   %B(abs(x-x_L) < xTol*max(1,abs(x)))=1;

   y   = ResultLP.y_k;
   v_k = ResultLP.v_k;
   cTx = c'*x;
   f_k = cTx;

% if ExitFlag==0 & abs(cTx-fDP)/max(1,abs(cTx)) > 1E-12 
%    disp('BIG DUALITY GAP');
%    format long
%    (cTx-fDP)/max(1,abs(cTx))
%    cTx
%    fDP
%    Iter
%    pause
% end

   if ExitFlag > 0 
      if PriLev >= 1
         disp('=+=+=+=+=+=+=+=')
         disp('=== cutplane:     Failure solving dual LP with simplex');
         disp('=+=+=+=+=+=+=+=')
      end

      ResultMIP.x_k  = x(1:n0);
% xprinte(x,'x:')
% ix=find(x < x_L)
% x(ix)-x_L(ix)
% ix=find(x > x_U)
% x_U(ix)-x(ix)
      ResultMIP.f_k  = f_k;
      ResultMIP.v_k  = v_k;
      ResultMIP.QP.B = B;
      ResultMIP.QP.y = y;
      ResultMIP.ExitFlag=ExitFlag;
      ResultMIP.ExitText=str2mat('Failure in dual LP solution. ',...
            ResultLP.ExitText);
      ResultMIP.QP.A   = A;
      ResultMIP.QP.b_U = b;
      ResultMIP.QP.b_L = b;
      ResultMIP.QP.x_L = x_L;
      ResultMIP.QP.x_U = x_U;
      ResultMIP.QP.c   = c;
      ResultMIP = endSolveMini(Prob,ResultMIP);
      return

   end

   Improve = cTx-cTx0;
   if PriLev > 1
      fprintf('cutplane: cTx =%40.20f Improvement: %15.7e\n',cTx,Improve)
   end
   cTx0=cTx;
   if PriLev > 2
      xprinte(x(1:n0),'x:');
   end
end % while

if PriLev >= 1
   fprintf('\n--- TOO MANY cutting plane ITERATIONS. ITER = %7.0f \n',Iter)
end
ResultMIP.x_k  = x;
ResultMIP.f_k  = f_k;
ResultMIP.v_k  = v_k;

ResultMIP.QP.A   = A;
ResultMIP.QP.b_U = b;
ResultMIP.QP.b_L = b;
ResultMIP.QP.x_L = x_L;
ResultMIP.QP.x_U = x_U;
ResultMIP.QP.c   = c;

ResultMIP.QP.B = B;
ResultMIP.QP.y = y;
ResultMIP.ExitFlag=1;
ResultMIP.ExitText=ExitText(1);
ResultMIP = endSolveMini(Prob,ResultMIP);

% ------------------------------------------------------------------------
function [B, Q, R, E, pRank, DelVar] = FullRank(B, A, Q, R, E, pRank,...
          pMax, epsRank) 
% ------------------------------------------------------------------------

nargin;

[m n]= size(A);

DelVar=[];
pRank0=pRank;

if pRank < pMax
   % Too few variables in basis or rank deficiency in linear constraints

   maxR=max(abs(diag(R)));
   [P,jE]=find(E);
   iN  = find(B == 1 | B==-1);

   if pRank==1
      iz=1:length(iN);
   else
      z=sum(abs(Q(:,1:pRank)'*A(:,iN)));
      [zz iz]=sort(z);
   end

   jMax=length(iz);

   j=0;

   while pRank < pMax & j < jMax
      j=j+1;
      k=iN(iz(j));
      BkOld=B(k);
      B(k)=0;
      iB = find(B==0);
      kB = find(k==iB);
      pRank1=pRank;

      [Q,R,P,pRank,maxR,absR] = InsertQR(Q, R, P, pRank, ...
           kB, A(:,k), epsRank);

      if pRank==pRank1
         % This column failed. Delete this column.

         B(k)=BkOld;
 
         [Q, R, P, pRank] = DeleteQR(Q, R, P, pRank, kB);

      else

         B(k)=0;

         if length(P) > m
            DelCol=P(m+1);
            [Q, R, P, pRank] = DeleteQR(Q, R, P, pRank, DelCol);
            k=iB(DelCol);
            % Put deleted variable on lower bound. Later check if OK.
            B(k)=1;
            % Save variable number for later check
            DelVar=[DelVar;k];

            %R=R(:,1:m);
            %DelCol=P(m+1:length(P));
            %while length(P) > m
            %   j=P(length(P));
            %   P(P > j)=P(P > j) - 1;
            %   P=P(1:length(P)-1);
            %end
            %k=iB(DelCol);
            %B(k)=1;
            %DelVar=[DelVar;k];
         
         end
      end
   end

if any(sort(P)~=(1:length(P))') 
   disp('lpSimplex (FullRank): Serious P error');
   disp('Bad input data, or some logical program error.');
   length(P)
   xprinti(P,'P:');
   pause
end

   if pRank > pRank0
      % Changes has been made to the QR decomposition
      E = sparse(P,1:m,ones(m,1),m,m);
      Q = Q(1:m,1:m);
      R = R(1:m,1:pRank);
   %else
   %   disp('NO CHANGES MADE IN FULLRANK')
   %   pause
   end

end

% ------------------------------
function Text = ExitText(Inform)
% ------------------------------

switch  Inform
   case 0
     Text = 'Optimal solution found';
   case 1
     Text = 'Maximal number of iterations reached';
   case 4
     Text = 'No solution found to LP relaxation';
end

% MODIFICATION LOG:
%
% 981114  hkh  Changed to new Result/Prob format for PhaseXSimplex, X=1,2
% 981115  hkh  Accessing Prob, not Prob1 and Prob2.
% 981119  hkh  Change field to Prob.QP.B
% 981124  hkh  Change dual LP call to DualSolve
% 990816  hkh  Revised completely for 2.0.
% 990913  hkh  Safeguard against nan and inf in x_0
% 040111  hkh  Add call to inisolve and endSolve, define optParam, field mLin
% 041123  hkh  Check MIP.VarWeight field to be able to handle LP problems
% 041222  med  Safeguard added for x_0
% 050421  hkh  Change lpSolve to lpSimplex, clean from DEBUG
% 051212  hkh  VarWht was wrongly used instead of VarWeight
% 051216  med  Removed print out
% 070222  hkh  Revise IntVars handling, use new format. Skip slacks as integers
% 070907  hkh  SolverLP/DLP picked from list, 1st with license, avoid GetSolver 
% 080606  med  Switched to iniSolveMini
% 080607  hkh  Switched to endSolveMini, use tomRun instead of tomSolve
% 080607  med  Switched to tomRunMini
