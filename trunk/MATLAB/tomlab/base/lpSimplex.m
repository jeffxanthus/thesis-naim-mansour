% lpSimplex.m
%
% Active set strategy (Simplex method) for Linear Programming.
%
% Minimization problem:
%
%
%        min    c' * x.  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%
% Equality equations: Set b_L(i)==b_U(i), for the ith equation
% Fixed    variables: Set x_L(i)==x_U(i), for the ith variable
%
% Solves Phase 1 problem (find feasible point) if c is empty
%
% function Result = lpSimplex(Prob)
%
% INPUT PARAMETERS
% Fields in Prob:
%   c:      The vector c in c'x. If isempty(c), only solve Phase1 problem,
%           finding a feasible point for the constraints
%   A:      The linear constraint matrix
%   b_L:    The lower bounds for the linear constraints
%   b_U:    The upper bounds for the linear constraints
%   x_L:    Lower bounds on x. If empty assumed to be 0.
%   x_U:    Upper bounds on x. If empty assumed to be Inf.
%           b_L, b_U, x_L, x_U must either be empty or of full length
%   x_0:    Starting point x
% PriLevOpt Printing level:
%           =0 No output; >0 Convergence results;
%           >1 Output every iteration of function value
%           >2 Output every iteration of change of basis
%           >3 Output every iteration of more information (scalar values)
%           >4 Output every iteration of more information including vectors
%  QP.B:    Active set B_0 at start.
%           2  = Fixed variables x(i)
%           1  = Include variable x(i) is in basic set.
%           0  = Variable x(i) is set on its lower bound
%           -1 = Variable x(i) is set on its upper bound
%           If EMPTY, lpSimplex finds the active set.
% -----------------------------------------------------------------------
% Warm Start basis variables:
% Prob.QP.UseHot   If > 0, Uses a warm start basis, either from file or from
%                    the structure Prob
% Prob.QP.HotFile If nonempty, a warm start basis are read from this file
%
% If Prob.QP.HotFile is empty, read warm start basis from structure Prob:
% x = Prob.QP.Hot.x; B = Prob.QP.Hot.B;
% Q = Prob.QP.Hot.Q; R = Prob.QP.Hot.R; E = Prob.QP.Hot.E;
%
% Prob.QP.HotFreq  If > 0, a basis are saved, on a file
% HotLPxxx.  xxx is a number (modulo counter).
%
% Prob.QP.HotN     The number of warm start basis files, xxx=mod(Freq,HotN),
%                    where xxx is a number from 0 to HotN-1
% -----------------------------------------------------------------------
%
% Prob.Solver.Alg  Rule to select new variables:
%           = 0 Minimum Reduced Cost, sort variables increasing. (DEFAULT)
%           = 1 Bland's Anti-cycling Rule
%           = 2 Minimum Reduced Cost, Dantzig's rule
%
% Fields used in Prob.optParam:
%
% wait      Wait flag, pause each iteration if set true
% MaxIter   Maximal number of iterations. max(10*dim(x),100) is DEFAULT.
% eps_f     Tolerance used to test convergence on the reduced costs
% eps_Rank  Rank tolerance used to test convergence on the reduced costs
% xTol      Tolerance to judge if x-values are close
% bTol      Feasibility tolerance for linear constraints
%
%
% OUTPUT PARAMETERS
% Structure Result. Fields used:
%   Iter     Number of iterations
%   ExitFlag Exit flag
%            == 0  => OK
%            == 1  => Maximal number of iterations reached. No bfs found.
%            == 2  => Unbounded feasible region.
%            == 3  => Rank problems (NOT USED)
%            == 4  => No feasible point found with Phase1
%            == 10 => Errors in input parameters
%   ExitTest Text string giving ExitFlag and Inform information
%   Inform   Empty except if ExitFlag=4: No of infeasible bounds/constraints
%   x_k      Solution
%   v_k      Lagrange parameters.
%            First n values corresponds to either lower or upper bounds.
%            The rest of the values corresponds to the linear constraints
%   QP.B     B  Optimal set. B(i)==0, include variable x(i) in basic set.
%            sum(B==0)==length(b_U)  holds. See QP.B as input.
%   p_dx     Search steps in x
%   alphaV   Step lengths for each search step
%   f_k      Function value c'*x
%   g_k      Gradient, =  c
%   Solver   lpSimplex
%   SolverAlgorithm  Description of method used
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written Nov 7, 1998.    Last modified Jun 6, 2008.

function ResultLP = lpSimplex(Prob)

if nargin < 1
   error('lpSimplex needs input structure Prob');
end

solvType=checkType('lp');

Prob=ProbCheck(Prob,'lpSimplex',solvType);

Prob = iniSolveMini(Prob);

% Define QP matrices and vectors from Prob structure

c      = Prob.QP.c(:);

[x, x_L, x_U, Prob, EMPTY, iFix, iFree] = GetX(Prob, length(c));

n = length(x_L);

[mA, Ax, bEqual, b_L, b_U, A] = LinearConstr(Prob);

bFree=isinf(b_L) & isinf(b_U);

m=mA;

optParam  = Prob.optParam;
PriLev    = Prob.PriLevOpt;       % Print level
IterPrint = optParam.IterPrint;   % Print short information each iteration

wait      = optParam.wait;        % Pause after printout if 1
epsRank   = optParam.eps_Rank;    % Rank test tolerance
MaxIter   = optParam.MaxIter;     % Maximal number of iterations
eps_f     = optParam.eps_f;       % Convergence tolerance on the reduced costs
xTol      = optParam.xTol;
bTol      = optParam.bTol;

B         = Prob.QP.B(:);

Ascale    = Prob.QP.Ascale;

% Scale the matrix A if Ascale is true
[Ascale, A, D1, D2, b_L, b_U, x, x_L, x_U, c] = ScaleA(A, Ascale, ...
    iFix, iFree, b_L, b_U, bTol, x, x_L, x_U, c);

if isempty(Prob.QP.UseHot),  Prob.QP.UseHot=0; end
if isempty(Prob.QP.HotFreq), Prob.QP.HotFreq=0; end
if isempty(Prob.QP.HotN),    Prob.QP.HotN=0; end

UseHot  = Prob.QP.UseHot;
HotFreq = Prob.QP.HotFreq;
HotN    = Prob.QP.HotN;

Freq=0;

if HotFreq > 0
   HotName='HotLP';
end

if length(b_U)~=m
   fprintf('A is %d by %d, but upper bounds b_L is %d\n',m,n,length(b_U));
   error('Illegal input to lpSimplex');
end

n      = max([length(c),size(A,2)]);
N      = mA+n;
lambda = zeros(N,1);

DEBUG=0;

ResultLP                 = ResultDef(Prob);
ResultLP.Prob            = Prob;
ResultLP.Solver          = 'lpSimplex';

[Alg,ResultLP.SolverAlgorithm] = SetAlg(Prob.Solver.Alg,'Simplex method.');
Prob.Solver.Alg=Alg;

if PriLev > 1
   fprintf('Run %s',ResultLP.SolverAlgorithm);
   fprintf('\n');
end

IndexRule = Prob.Solver.Alg;

ResultLP.f_0 = [];

n0     = n;
if issparse(A)
   A = [A -speye(m,m)];
else
   A = [A -eye(m,m)];
end
n = n+m;

% Restrict output
if n < 40
   mP=n;
else
   mP=40;
end

%rTol=[xTol*ones(n0,1);bTol*ones(m,1)];
rTol=[bTol*ones(n0,1);bTol*ones(m,1)];

if isempty(c)
   c     = zeros(n,1);
   SolveFP=1;
   Phase = 1;
else
   c     = [c;zeros(m,1)];
   SolveFP=0;
   Phase = 2;
end
c0    = c;
bl = [x_L;b_L];
bu = [x_U;b_U];

% Detect bound infeasible problem
if any( bl-bu > rTol .* max(1,abs(bl+bu)))
   ResultLP.ExitFlag=4;
   ResultLP.ExitText=ExitText(4);
   ResultLP.Inform=sum( bl-bu > rTol .* max(1,abs(bl+bu)));
   ResultLP.QP.B=[];
   if PriLev > 1
      disp('lpSimplex detected infeasible problem');
      [bl bu]
      if wait, pause; end
   end
   ResultLP = endSolve(Prob,ResultLP);
   return;
end

if ~EMPTY & ~UseHot
   % Only accept x with length n (n+m) or n0
   if ~(length(x)==n | length(x)==n0)
      EMPTY=1;
   end
   if any(isinf(x))
      % Do not accept inf values
      EMPTY=1;
   end
end

if ~isempty(B) & ~UseHot
   % Use given B
   if ~isempty(iFix)
      B(iFix)=2;
   end
   m0=sum(B==0);
   if m0==0 | m0 > m % User error. Skip using B
      B=[];
   end
end

if UseHot
   [x,B,Q,R,E,maxR,pRank] = HotRead(Prob,epsRank,PriLev);
else
if mA==0
   x
   %x=min(max(0,x_L),x_U);
   %x
   PriLev=5 
   wait
end
   [x, B, Q, R, E, pRank, Stop] = GetXB(x, B, bl, bu, iFree, bFree,...
             bEqual, iFix, A, b_L, b_U, m, n0, EMPTY, rTol, epsRank, Prob);
   if Stop==1
      B = ~((x > x_L + xTol) & (x < x_U - xTol) & (x_L~=x_U));
      B = double(B);
      B(find(x >= x_U-xTol))=-1;
      ResultLP = SaveResult(2000, Ascale, x, zeros(m,1), D1, D2, B, Q, R, E,...
          1, n0, m, 0, [], [], pRank, x_L, x_U, ResultLP);
      ResultLP = endSolve(Prob,ResultLP);
      return;
   elseif Stop==2
      if PriLev > 0
         disp('lpSimplex: Illegal initial logical basis B given as input')
         xprinti(B(1:mP),'B:');
         xprinte(xScale(x, Ascale, D2, mP),'x:');
      end
      ResultLP.x_k=xScale(x, Ascale, D2, n0);
      ResultLP.QP.B=B(1:n0);
      ResultLP.ExitFlag=10;
      ResultLP.ExitText=ExitText(10);
      ResultLP = endSolve(Prob,ResultLP);
      return
   elseif Stop==3
      disp('Illegal initial x to lpSimplex')
      xprinte(xScale(x, Ascale, D2, mP),'x:');
      ResultLP.x_k=xScale(x, Ascale, D2, n0);
      ResultLP.QP.B=B(1:n0);
      ResultLP.ExitFlag=11;
      ResultLP.ExitText=ExitText(11);
      ResultLP = endSolve(Prob,ResultLP);
      return
   end
   
end

[iB,N_L,N_U,iN] = Bsets(B);
%iB  = find(B==0);
%iN  = find(B == 1 | B==-1);


if pRank < m
   % Generate a full rank basis, if possible
   pRank0=pRank;
   [B, Q, R, E, pRank, DelVar] = FullRank(B, A, Q, R, E, pRank, m, epsRank);
   if PriLev > 1 & pRank0 < pRank
      fprintf('FullRank: Increased rank from %d to %d.',pRank0,pRank);
      fprintf(' Rows m = %d.\n',m);
   end
   if ~isempty(DelVar)
      x(DelVar)=bl(DelVar);
      B(isinf(x))=-1;
      x(isinf(x))=bu(isinf(x));
      x(isinf(x))=0;            % Free variables set to 0
   end
   if pRank < m
      % Correct the values for the variables in the base which will not
      % be moved during the iterations
      [P,jE]=find(E);
      iB = Bsets(B);
      ix=iB(P(pRank+1:m));
      x(ix)=bl(ix);             % First set all nonmoving variables to zero
      x(isinf(x))=bu(isinf(x)); % Set to upper bound if lower is inf
      x(isinf(x))=0;            % Free variables set to 0
   end

   [iB,N_L,N_U,iN] = Bsets(B);
   %iB  = find(B==0);
   %iN  = find(B == 1 | B==-1);
end

if length(iB) > m
   % Remove variables from basis, based on the result from the QR
   [P,jE]=find(E);
   % Remove extra unnecessary columns
   for i=length(iB):-1:m+1
       % Drop number i. Renumber P
       j = P(i); 
       k = iB(j);
       if abs(x(k)-x_U(k)) < xTol
          x(k)=x_U(k);
          B(k)=-1;
       else
          x(k)=x_L(k);
          B(k)=1;
       end
       P(P > j)=P(P > j) - 1;
       iB = Bsets(B);
   end
   P=P(1:m);
   R=R(:,1:m);
   iN  = find(B == 1 | B==-1);
   E=sparse(P(:),[1:m]',ones(m,1),m,m);
end

x_0 = x;  % Save x_0. Needed if rank problems in base

if any(isinf(x)) 
   PriLev=3
   wait=1
   if PriLev > 0
      [bl x bu]
      disp('inf at start')
      if wait, pause; end
   end
end

%PriLev=5
%wait =1

if pRank < m
   % Correct the values for the variables in the base which will not
   % be moved during the iterations. They should be on some bound,
   % which may not be 0 always.
   [P,jE]=find(E);
   iB = Bsets(B);
   ix=iB(P(pRank+1:m));
   Corr=A(:,ix)*x_0(ix);
else
   Corr=zeros(m,1);
end
if isempty(iFix)
   x(iB) = tomsol(6, Q, R, E, pRank, -A(:,iN)*x(iN)-Corr);
else
   if isempty(iN)
      x(iB) = tomsol(6, Q, R, E, pRank, -A(:,iFix)*x(iFix)-Corr);
   else
      x(iB) = tomsol(6, Q, R, E, pRank, -A(:,iN)*x(iN)-A(:,iFix)*x(iFix)-Corr);
   end
end
if pRank < m
   % Put the nonchanging values to correct values.
   x(ix)=x_0(ix); 
end

ResultLP.x_0 = x(1:n0);
ResultLP.f_0 = c(1:n0)'*x(1:n0);
ResultLP.g_k = c(1:n0);

[pOld,f_kOld,rSumOld,Cycle,doQR,qrFull,sinceQR,ixOld,ixL,ixU,y ...
  Pen,pRank0,rAx,N_L,N_U,Iter] = InitVars(c,bl,bu,rTol,m,n0,B,pRank,PriLev);

while Iter < MaxIter

   Iter = Iter + 1;
   ResultLP.Iter     = Iter;

   if Iter > 1 & Phase==2
      % Compute slack variable values
      Ax = A(:,1:n0)*x(1:n0);
      rS=x(n0+1:n0+m) - Ax;
%     xprinte(rS,'rS');
      rSnorm=norm(rS);
      maxrS=max(rS);
      % Check the feasibility
      if maxrS > 1E-10 * max(1,abs(Ax))
         doQR=1;

         format compact
         rSnorm
         maxrS
         Iter
%if maxrS > 1E-3
%   PriLev=5
%   wait=1
%   pause
%end
      end
   end

   if PriLev > 2
      xprinte(xScale(x, Ascale, D2, mP),'x:');
   end

   if any(isinf(x))
      disp('lpSimplex: Some x is inf')
      pause
      %xprint(x,'x:');
      %Iter
      %PriLev=3
      %wait=1
      %pause
   end
   if pRank < pRank0 & PriLev > 1
      fprintf('Rank has decreased. Iter %d. Rank %d. ',Iter,pRank);
      fprintf('Rank at start %d. m = %d\n',pRank0,m);
   end


   % Check if problem is feasible
 
   %xprinte(x(1:10),'x0:');
   %r =max(0,(bl-x)./max(1,abs(x))) + max(0,(x-bu)./max(1,abs(x)));
   %plot(r(r>0))
   %pause

   %r =max(0,(bl-x)./max(1,abs(x))-rTol) + max(0,(x-bu)./max(1,abs(x))-rTol);
   rL =max(0,bl-x);
   rU =max(0,x-bu);
   r =max(0,rL+rU);
   rSum = sum(r);
   rx = rTol(1:n0).*max(1,abs(x(1:n0)));
   ix = find(abs(x(1:n0)) <= xTol | x(1:n0) < x_L(1:n0));
   if ~isempty(ix)
      rx(ix) = mean(rx);
   end
   
   if any(r(1:n0) > rx) | any(r(n0+1:n) > rAx)

if PriLev > 4
xprinte(r,'r:');
end
if PriLev > 1 
idx=find(r>eps);
xprinte(r(idx(1:min(6,length(idx)))) ,'r:');
%pause
end

      % Do Phase I to restore feasibility
      Phase=1;
      ixL=rL > rTol .* max(1,abs(x));
      ixU=rU > rTol .* max(1,abs(x));
      P1Var=find(ixL & B ~= 0);
      if sum(ixU) > 0
         P1Var=[P1Var;-find(ixU & B ~= 0)];
      end

      P1VarBug(P1Var,bl,bu,x,B,iB,Iter,n0,m,pRank,PriLev);
      
      c = c0 + Pen * [ ixU - ixL];

      %bu(ixL & isinf(bu))=bl(ixL & isinf(bu));

      bl(ixL)=-Inf;

      %bl(ixU & isinf(bl))=bu(ixU & isinf(bl));

      bu(ixU)=Inf;
   
   elseif Phase==1
      if SolveFP
         % Only Phase 1 problem.
         ResultLP = SaveResult(1000, Ascale, x, y, D1, D2, B, Q, R, E,...
             Iter, n0, m, 0, [], N_U, pRank, x_L, x_U, ResultLP);

         if PriLev > 0
            fprintf('\n=======******-------------------------------------')
            fprintf('-----------------------\n')
            fprintf('lpSimplex: Phase 1 converged in %d iterations\n',Iter)
            fprintf('=======******-------------------------------------')
            fprintf('-----------------------\n')
            if PriLev > 1
               fprintf('B: %d basic variables:\n',length(iB));
               xprinti(iB,'B:  ');
               if ~isempty(iFix)
                  fprintf('F: %d fixed variables:\n',length(iFix));
                  xprinti(iFix,'F:  ');
               end
               if ~isempty(N_L)
                  fprintf('N_L: %d non-basic',length(N_L));
                  fprintf(' lower bound variables:\n');
                  xprinti(N_L,'N_L:');
               end
               if ~isempty(N_U)
                  fprintf('N_U: %d non-basic',length(N_U));
                  fprintf(' upper bound variables:\n');
                  xprinti(N_U,'N_U:');
               end
            end
         end
         if PriLev > 2
            xprinte(xScale(x, Ascale, D2, mP),'x:');
            disp('Base-B Lower-Bound  [x;slack s]  Upper-Bound');
            [B bl x bu]
         end
         if wait & PriLev > 1
            disp('Press any key to continue ...')
            pause
         end
         ResultLP = endSolve(Prob,ResultLP);
         return
      end
      % Reset lower bound, upper bound, and objective function
      bl = [x_L;b_L];
      bu = [x_U;b_U];
      [iB,N_L,N_U,iN] = Bsets(B);
      c  = c0;
      Cycle=0;
      pOld=0;
      f_kOld  = NaN;
      rSum    = 0;
      rSumOld = 0;
      ixL    = [];
      ixU    = [];
      Phase=2;
   end

   %===== Step 2:  Compute the objective function value c^T * x

   f_k = c0' * x;

   % ===== Step 2 to 8 in the simplex algorithm ===================

   % B(j) = 0 if variable j in basis B.  Else B(j) = 1 (low) or -1 (upp).
   % IndexRule == 0 Minimum Reduced Cost
   %           == 1 Blands rule
   %           == 2 Minimum Reduced Cost, Dantzig algorithm
   % p              Variable to include in basis
   % q              Variable to exclude in basis
   % alpha          Step length
   % d              Search direction
   % f_k            c'*x, objective function
   % y              Dual variables,Lagrangian multipliers to linear constraints
   % cHat           Reduced costs, Lagrangian multipliers to bounded variables



   if f_k == f_kOld & rSum==rSumOld
      Cycle=Cycle + 1;
   else
      Cycle=0;
   end
   f_kOld  = f_k;
   rSumOld = rSum;

   if PriLev > 1
      fprintf('\n--------------------------------------------------')
      fprintf('-----------------------\n')
      fprintf('Phase %d.',Phase)
      fprintf(' Iter: %5.0f. ',Iter)
      fprintf('f_k: %25.16f',f_k)
      if rSum > 0
         fprintf(' Infeas: %10.5e',rSum)
      end
      fprintf('\n')
      fprintf('--------------------------------------------------')
      fprintf('-----------------------\n')
   elseif IterPrint
      fprintf('Phase %d.',Phase)
      fprintf(' Iter: %5.0f. ',Iter)
      fprintf('f_k: %25.16f',f_k)
      if rSum > 0
         fprintf(' Infeas: %10.5e',rSum)
      end
      fprintf('\n')
   end
   if PriLev > 2
      fprintf('B: %d basic variables:\n',length(iB));
      xprinti(iB,'B:  ');
      if ~isempty(iFix)
         fprintf('F: %d fixed variables:\n',length(iFix));
         xprinti(iFix,'F:  ');
      end
      if ~isempty(N_L)
         fprintf('N_L: %d non-basic lower bound variables:\n',length(N_L));
         xprinti(N_L,'N_L:');
      end
      if ~isempty(N_U)
         fprintf('N_U: %d non-basic upper bound variables:\n',length(N_U));
         xprinti(N_U,'N_U:');
      end
   end
   if PriLev > 4
      xprinte(xScale(x, Ascale, D2, mP),'x:');
   end

   %===== Step 3:  Compute shadow prices =========================
   % y = A(:,iB)' \ c(iB);
   % QR - factorization already computed: [Q R E] = qr(full(A(:,iB)));

   y = tomsol(7,Q, R, E, pRank, c(iB(1:m)));

   if PriLev > 4
      fprintf('Shadow prices:\n');
      xprinte(y,'y: ');
   end

   %===== Step 4:  Compute reduced costs =========================
   % Only compute reduced costs for nonbasic variables
   %cHat = zeros(n,1);
   %cHat(iN) = c(iN) - A(:,iN)' * y;

   % Compute reduced costs for all variables. 
   cHat = c - A' * y;
   if ~isempty(N_U)
      cHat(N_U) = -cHat(N_U); % Change signs on upper bound reduced costs
   end
   % Check that reduced costs for basic variables sufficiently zero
   if max(abs(cHat(iB))) > bTol
      % Too bad numerical accuracy. Recompute QR decomposition
      if PriLev > 1
         disp('Bad numerical accuracy. Recompute QR')
      end
      %plot(cHat(iB))
      %xprinte(cHat(iB),'cHat:');
      %pause
      doQR=1;
   end

   cH=cHat;
   % Remove fixed variables from min cost test
   cH(find(~isinf(bl) & ~isinf(bu) & (abs(x-bl) < rTol .* max(1,abs(x))) & ...
            (abs(bl-bu) < rTol .* max(1,abs(x)))))=0;

   % Set reduced cost 0 for variables in basis
   cH(iB) = 0;

   if PriLev > 3
      if ~isempty(N_L)
         fprintf('Reduced costs for lower bounds:\n');
         xprinte(cHat(N_L),'  ');
      end
      if ~isempty(N_U)
         fprintf('Reduced costs for upper bounds:\n');
         xprinte(cHat(N_U),'  ');
      end
   end

   %===== Step 5:  Check if minimum, Reduced costs all nonnegative =========

   [minCost imC] = sort(cH);

   qVar=0;
   if minCost(1) > -eps_f & Phase==1
      if Pen > 1E12 | minCost(1) >= -eps

         if PriLev > 0
            % Convergence of Phase 1 iteration WITHOUT Feasibility
            disp('=======================================');
            disp('lpSimplex: NO FEASIBLE SOLUTION FOUND!!!!');
            disp('=======================================');
            fprintf(' Infeas: %10.5e\n',rSum)
         end
         if PriLev > 2
            disp('Lower-bound-x x Upper-bound-x')
            [x_L x(1:n0) x_U]
            if n0 > 20 & wait
               pause
            end
            disp('Lower-Linear slack-s A*x Upper-Linear-bound')
            [b_L x(n0+1:n) A*x b_U]
            pause
         end
         if wait & PriLev > 1
            disp('Press any key to continue ...')
            pause
         end
         ResultLP = SaveResult(4, Ascale, x, y, D1, D2, B, Q, R, E,...
             Iter, n0, m, 0, cHat, N_U, pRank, x_L, x_U, ResultLP);

         % Return the residual vector + logical vector
         ResultLP.r_k=[r(1:n0) > rx ; r(n0+1:n) > rAx];
         ResultLP.QP.Hot.x = x;
         ResultLP.QP.Hot.B = B;
         ResultLP.QP.Hot.Q = Q;
         ResultLP.QP.Hot.R = R;
         ResultLP.QP.Hot.E = E;

         ResultLP = endSolve(Prob,ResultLP);
         return;
      else   % Increase penalty
         Pen=Pen*10;
         minCost=minCost*1E15;
      end
   end
   if minCost(1) > -eps_f
      if PriLev >= 1
         fprintf('\n=======******-------------------------------------')
         fprintf('-----------------------\n')
         fprintf('Minimum Phase %d. ',Phase)
         fprintf(' Iteration: %5.0f. ',Iter)
         fprintf('f_k: %30.20f\n',f_k)
         fprintf('\n=======******-------------------------------------')
         fprintf('-----------------------\n')
      end
      if wait & PriLev > 1
         disp('Press any key to continue ...')
         pause
      end
      if PriLev >= 3
         fprintf('\n')
         xprinte(xScale(x, Ascale, D2, mP),'x:');
      end
      ResultLP = SaveResult(0, Ascale, x, y, D1, D2, B, Q, R, E,...
          Iter, n0, m, f_k, cHat, N_U, pRank, x_L, x_U, ResultLP);
      ResultLP = endSolve(Prob,ResultLP);
      return;
   end

   %===== Step 6:  Choose index for variable to include in the new basis =====
   useBland=1;
   if Cycle > n
      if PriLev >= 0
         fprintf('lpSimplex: Cycling! Iteration %d\n',Iter)
      end
      if IndexRule==1 | IndexRule==3
         useBland=0;
         if PriLev > 0
            disp('lpSimplex: Bland''s rule does not work!')
         end
      else
         % Switch to Bland's rule
         IndexRule=1;
      end
   end
   nCand = sum(minCost <= -eps_f);
   if (IndexRule == 1 | IndexRule == 3) & useBland    		
      % Choose index using Blands regel

      [q iq] = sort(imC(1:nCand));

      if pOld==q(1) & nCand > 1
         q = q(2);
         mc= minCost(2);
      else
         q = q(1);
         mc= minCost(1);
      end
   else                   		
      % Choose index using Minimal Reduced Cost
      if pOld==imC(1) & nCand > 1
         q = imC(2);
         mc= minCost(2);
      else
         q = imC(1);
         mc= minCost(1);
      end
   end

   if isempty(N_L)
      LOW = 0;
   else
      LOW  = any(q==N_L);
   end
   if LOW
      qHat= find(q==N_L);
   else
      qHat= find(q==N_U);
   end

   if PriLev > 2
      SS='UL';
      fprintf('x_q from N_%s',SS(LOW+1));
      fprintf(' added to basic set B:         ');
      fprintf('q = %5.0f, q^ = %5.0f',q,qHat)
      fprintf('.  minCost = %20.15f\n',mc);
   end

   %===== Step 7: Compute the search direction d =========================
   % d = A(:,iB) \ (-+A(:, q)); if q is in Lower och Upper bound set

   if LOW
      d = tomsol(6, Q, R, E, pRank, -A(:,q));
      dErr=sum(abs(A(:,iB)*d+A(:,q)));
   else
      d = tomsol(6, Q, R, E, pRank, A(:,q));
      dErr=sum(abs(A(:,iB)*d-A(:,q)));
   end

   %if PriLev > 4
   %   xprinte(d,'d:');
   %end

   if dErr > 1E-8
      if PriLev > 1
         disp('Residual for search direction too large. Force QR rebuild')
         fprintf('|A(:,iB)*d +- A(:,q)| = %23.10e\n',dErr);
         if ~isempty(iFix)
            xprinte(A(:,iFix)*x(iFix),'AxFix');
         end
         fprintf('Feasibility? 0 = |A*x - Is | = %23.10e\n',norm(A*x));
      end
      doQR=1;
   end

   %===== Step 8: Compute steplength alpha 
   %              and variable p to exclude from basis
   [alpha, p, pHat, LOW, Implicit] = alphaStep(m, iB, x, d, bl, bu, ...
    x_L, x_U, b_L, b_U, ixL, ixU, LOW, IndexRule, Phase, rTol, q, qVar, PriLev);

   if Phase==1
      bl = [x_L;b_L];
      bu = [x_U;b_U];
   end

   if isinf(alpha) | alpha > 1E10
      %  if alpha > 1E10  --- Assume unbounded for such big steps
      ExitFlag = 2;
      % Problem has unbounded feasible region!

      % Check angles

      %w=A(:,iB)'*A(:,q)/norm(A(:,q));
      %for j=1:length(w), w(j)=w(j)/norm(A(:,iB(j))); end
      %w=180*acos(abs(w))/pi;
      %for j=1:length(w)
      %    if w(j) < 1
      %       fprintf('WARNING! ');
      %       fprintf('Near linear dependence between var q=%d',q);
      %       fprintf(' and var p=%d. Angle %7.3f\n',iB(j),w(j));
      %    end
      %end
      ResultLP = SaveResult(2, Ascale, x, y, D1, D2, B, Q, R, E,...
          Iter, n0, m, f_k, cHat, N_U, pRank, x_L, x_U, ResultLP);
      ResultLP.ExitText=ExitText(2);
      ResultLP = endSolve(Prob,ResultLP);
      return;
   end

   % ===== Step 9 Update of x and B =========================


   x(iB) = x(iB) + alpha * d;
   if abs(x(p)-bl(p)) > bTol*max(1,abs(x(p))) & ...
      abs(x(p)-bu(p)) > bTol*max(1,abs(x(p)))
      disp('Illegal step. x(p) should be on some bound:')
      fprintf('x(p) %20.15f x_L %20.15f x_U %20.15f\n',x(p),bl(p),bu(p));
   end

   xqOld = x(q);
   %x = xUpdate(x, q, bl, bu, alpha, Implicit, qVar, LOW);
   if alpha > 0 & ~Implicit
      if LOW
         x(q) = bl(q) + alpha;
      else
         x(q) = bu(q) - alpha;
      end
   elseif alpha > 0 & qVar~=0
      if LOW
         x(q) = x(q) + alpha;
      else
         x(q) = x(q) - alpha;
      end
   end

   [iB, N_L, N_U, B, x, doQR] = UpdateB(p, pHat, q, qHat, rTol, ...
    iB, N_L, N_U, B, x, bl, bu, IndexRule, LOW, doQR, PriLev) ;

if length(iB) > m
   disp('lpSimplex: Error! Number of variables in basis are too large');
   xprinti([m;Iter;iB(:)],'m,Iter,iB')
   pause
end
   iN = find(B==1 | B==-1);

   pOld=p;  % Save the last variable that leaves the basis, to avoid it
            % coming back directly, if possible.

   if pHat ~= 0
      if (mod(Iter,qrFull)==0 & qrFull > 0) | (doQR > 0 & sinceQR > 10)  
         if doQR > 0
            if PriLev > 1
               fprintf('FORCED QR!!!\n');
            end
            doQR=0;
         else
            if PriLev > 1
               fprintf('Rebuild QR!!! Iter %d qrFull %d\n',Iter,qrFull);
            end
         end
         % Make QR factorization

         [Q, R, E, pRank] = ComputeQR(A(:,iB), epsRank);

         sinceQR=0;
      else
         sinceQR=sinceQR+1;

         [P,jE]=find(E);

         % Delete column pHat, variable p, from QR-decomposition
 
         [Q, R, P, pRank] = DeleteQR(Q, R, P, pRank, pHat);

         % Insert column for variable q, into QR-decomposition
         % k is new column position in iB.

         k=find(q==iB);
         if isempty(k)
            error('lpSimplex: Fatal Error. Empty variable');
         end

         [Q,R,P,pRank,maxR,absR] = InsertQR(Q, R, P, pRank, k, A(:,q), epsRank);

if any(sort(P)~=[1:m]')
   disp('lpSimplex: Serious P error');
   disp('Bad input data, or some logical program error.');
   xprinti(P,'P:');
   pause
end
         E=sparse(P(:),[1:m]',ones(m,1),m,m);
         if pRank < pRank0
            % Try to increase rank up to initial rank
            [B, Q, R, E, pRank, DelVar] = FullRank(B, A, Q, R, E, pRank,...
             pRank0, epsRank);
            if ~isempty(DelVar)
               x(DelVar)=bl(DelVar);
               B(isinf(x))=-1;
               x(isinf(x))=bu(isinf(x));
            end
            [iB,N_L,N_U,iN] = Bsets(B);
            %disp('RECOMPUTE !!!')
            %PriLev=2;
            vv= -A(:,iN)*x(iN)-Corr;
            if isempty(iFix)
               x(iB) = tomsol(6, Q, R, E, pRank, vv);
               %x(iB) = tomsol(6, Q, R, E, pRank, -A(:,iN)*x(iN)-Corr);
            else
               vv=vv-A(:,iFix)*x(iFix);
               x(iB) = tomsol(6, Q, R, E, pRank,vv);
               %x(iB) = tomsol(6, Q, R, E, pRank, ...
               %               -A(:,iN)*x(iN)-A(:,iFix)*x(iFix)-Corr);
            end
         end
      end
   end % if pHat ~= 0

   if wait & PriLev > 1
      disp('Press any key to continue ...')
      pause
   end
   if mod(Iter,HotFreq)==0
      Freq=Freq+1;
      %save([HotName num2str(mod(Freq,HotN))],'x','B','Q','R','E');
      save([HotName num2str(mod(Freq,HotN))],'x','B');
   end
end

if mod(Iter,HotFreq)~=0 & HotFreq > 0
   Freq=Freq+1;
   save([HotName num2str(mod(Freq,HotN))],'x','B','Q','R','E');
end
if PriLev >= 0
   fprintf('\nlpSimplex: TOO MANY ITERATIONS --- ITER =%7.0f\n',Iter)
   fprintf('\nEither cycling or some numerical problems occurred\n')
   %fprintf('\nPress ENTER to continue\n')
   %pause
end

ResultLP = SaveResult(1, Ascale, x, y, D1, D2, B, Q, R, E,...
   Iter, n0, m, f_k, cHat, N_U, pRank, x_L, x_U, ResultLP);
% Also return the residual vector, in case of Phase 1
ResultLP.r_k  = r;
ResultLP = endSolve(Prob,ResultLP);

% ------------------------------------------------------------------------

function [iB, N_L, N_U, B, x, doQR] = UpdateB(p, pHat, q, qHat, rTol, ...
          iB, N_L, N_U, B, x, x_L, x_U, IndexRule, LOW, doQR, PriLev) 
% ------------------------------------------------------------------------

nargin;

B = double(B);

B(q) = 0;
if abs(x(p) - x_L(p)) <= rTol(p) * max(1,abs(x(p)))
   SETLOW=1;
elseif abs(x(p) - x_U(p)) <= rTol(p) * max(1,abs(x(p)))
   SETLOW=0;
else
   % Bad accuracy. Try to do new QR
   doQR=1;
   if abs(x(p) - x_L(p)) <= abs(x_U(p)-x(p))
      SETLOW=1;
   else
      SETLOW=0;
   end
   if PriLev >= 0
      disp('lpSimplex: Variable is NOT close to bound')
      xprinte([x_L(p) x(p) x_U(p) rTol(p)],'Trouble:')
[x_L x x_U]
      pause
   end
end

if SETLOW
   x(p) = x_L(p); % Set identically to bound, avoid round off difficulties 
   B(p) = 1;
   if IndexRule ==0
      [iB,N_L,N_U] = Bsets(B);
   elseif IndexRule ==1
      [iB,N_L,N_U] = Bsets(B);
   elseif IndexRule ==2
      N_L(qHat) = p;
      if pHat > 0, iB(pHat)  = q; end
   elseif IndexRule ==3
      if pHat > 0, iB  = [iB(1:pHat-1);iB(pHat+1:length(iB));q]; end
      if LOW
         N_L = [N_L(1:qHat-1);N_L(qHat+1:length(N_L));p];
      else
         N_U = [N_U(1:qHat-1);N_U(qHat+1:length(N_U))];
         N_L = [N_L;p];
      end
   end
else
   x(p) = x_U(p); % Set identically to bound, avoid round off difficulties 
   B(p) = -1;
   if IndexRule ==0
      [iB,N_L,N_U] = Bsets(B);
   elseif IndexRule ==1
      [iB,N_L,N_U] = Bsets(B);
   elseif IndexRule ==2
      N_U(qHat) = p;
      if pHat > 0, iB(pHat)  = q; end
   elseif IndexRule ==3
      if pHat > 0, iB  = [iB(1:pHat-1);iB(pHat+1:length(iB));q]; end
      if LOW
         N_L = [N_L(1:qHat-1);N_L(qHat+1:length(N_L))];
         N_U = [N_U;p];
      else
         N_U = [N_U(1:qHat-1);N_U(qHat+1:length(N_U));p];
      end
   end
end

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

if any(sort(P)~=[1:length(P)]') 
   disp('lpSimplex (FullRank): Serious P error');
   disp('Bad input data, or some logical program error.');
   %size(R)
   %size(E)
   m
   length(P)
   xprinti(P,'P:');
   pause
end


   if pRank > pRank0
      % Changes has been made to the QR decomposition
      E = sparse(P,[1:pRank],ones(pRank,1),pRank,pRank);
      Q = Q(1:m,1:m);
      R = R(1:m,1:pRank);
   %else
   %   disp('NO CHANGES MADE IN FULLRANK')
   %   pause
   end

   % Logical handling fix
   B = double(B);
   
end

% ------------------------------------------------------------------------
function [alpha, p, pHat, LOW, Implicit] = alphaStep(m, iB, x, d, bl, bu,...
   x_L, x_U, b_L, b_U, ixL, ixU, LOW, IndexRule, Phase, rTol, q, qVar, PriLev) 
% ------------------------------------------------------------------------

%===== Step 8: Compute steplength alpha and variable p to exclude from basis

nargin;

Step = Inf*ones(m,1);
L=find(~isinf(bl(iB)) & d < -rTol(iB) .* max(1,abs(x(iB))));
U=find(~isinf(bu(iB)) & d >  rTol(iB) .* max(1,abs(x(iB))));
if ~isempty(L)
   Step(L)=(bl(iB(L))-x(iB(L)))./d(L);
end
if ~isempty(U)
   Step(U)=(bu(iB(U))-x(iB(U)))./d(U);
end
alpha=min(Step);

%if Phase==1 & (alpha==inf | alpha==0) & 0
%   [bl x bu]
%   [d Step]
%   alpha
%   pause
%end
 
if ~isinf(alpha)
   ip = find(alpha==Step);
   if IndexRule == 0 | IndexRule==2
      % In case of ties, select the column with the highest index 
      % to leave 
      pHat=ip(length(ip));
   elseif IndexRule==1 | IndexRule == 3  
      % In case of ties, select the column with the lowest index 
      % to leave (Blands rule)
      pHat=ip(1);
   end
end
if ~isinf(bu(q)) & ~isinf(bl(q)) 
   dist=bu(q)-bl(q);
   if dist <= alpha & dist > rTol(q) * max([1,abs(bu(q)),abs(bl(q))])
      % The released variable is going directly to its other bound
      % In case of ties, use this variable first
      alpha = dist;

      % disp('The released variable goes directly to its other bound');
      % alpha
      pHat = 0;
   end
end
if isinf(alpha)
   % Check for degeneracy
   il=~isinf(bl(iB)) & ( abs(x(iB)-bl(iB)) <= ...
       rTol(iB) .* max(1,abs(bl(iB))) )  & ...
      (abs(d) <= rTol(iB) .* max(1,abs(bl(iB))));
   iu=~isinf(bu(iB)) & ( abs(x(iB)-bu(iB)) <= ...
       rTol(iB) .* max(1,abs(bu(iB))) )  & ...
   (abs(d) <= rTol(iB) .* max(1,abs(bu(iB))));
   if any(il | iu)
      pHatx=find(il);
      if ~isempty(pHatx)
         pHat=pHatx(1);
      else
         pHatx=find(iu);
         pHat=pHatx(1);
      end
      alpha=0;
   else
      alpha=Inf;
   end
end
if Phase==1
   bl = [x_L;b_L];
   bu = [x_U;b_U];
end
if Phase == 1
   Step = Inf*ones(m,1);
   L=find(ixL(iB) & d >  rTol(iB) .* max(1,abs(x(iB))));
   U=find(ixU(iB) & d < -rTol(iB) .* max(1,abs(x(iB))));
   if ~isempty(L)
      Step(L)=(bl(iB(L))-x(iB(L)))./d(L);
   end
   if ~isempty(U)
      Step(U)=(bu(iB(U))-x(iB(U)))./d(U);
   end
   alpha0=min(Step);
   if alpha0 < alpha
      alpha=alpha0;
      pHat=find(alpha==Step);
      pHat=pHat(1);
     

   end

end
Implicit=0;
if isinf(alpha) & Phase == 1
   % Use implicit bounds
   if qVar > 0
      alpha=bl(qVar)-x(qVar);
      LOW=1;
      pHat=0;
      Implicit=1;
   elseif qVar < 0
      LOW=0;
      qVar=abs(qVar);
      alpha=x(qVar)-bu(qVar);
      pHat=0;
      Implicit=1;
   end
end

if isinf(alpha) | alpha > 1E10
   %  if alpha > 1E10  --- Assume unbounded for such big steps
   % Problem has unbounded feasible region!
   if PriLev >= 1
      fprintf('\n\nlpSimplex: ');
      fprintf('Phase %d simplex: Unbounded feasible region\n',Phase);
      if any((x(iB) - bl(iB) <=  rTol(iB) .* max(1,abs(x(iB))) | ...
              x(iB) - bu(iB) >= -rTol(iB) .* max(1,abs(x(iB)))) ...
              & (abs(d) <= rTol(iB) .* max(1,abs(x(iB)))))
         fprintf('Phase %d simplex: Problem is degenerated!\n',Phase);
      end
   end
   p=[];
   pHat=[];
else
   if pHat > 0 
      p = iB(pHat); 
   else 
      p = q;
   end

   if PriLev > 2
      fprintf('Step to nearest constraint is alpha = %20.10f \n',alpha)
      fprintf('x_p deleted from basic set B:              ')
      fprintf('p = %5.0f, p^ = %5.0f',p,pHat)
      fprintf('. x(p)= %20.15f\n',x(p));
   end
end

% ------------------------------------------------------------------------
function [x, x_L, x_U, Prob, EMPTY, iFix, iFree] = GetX(Prob, n)
% ------------------------------------------------------------------------

x = Prob.x_0(:);  
if any(isnan(x) | isinf(x)), x=[]; end
if isempty(x)
   EMPTY=1;
elseif all(x==0)
   % Treat x zero as empty solution
   EMPTY=1;
else
   EMPTY=0;
end

x_L = Prob.x_L(:);
x_U = Prob.x_U(:);

n = max([length(x_L),length(x_U),n,length(x)]);

if isempty(x_U),x_U= Inf*ones(n,1); end
if isempty(x_L),x_L=zeros(n,1); end

if length(x_L) < n
   if length(x_L)==1 
      x_L=x_L*ones(n,1);
   else 
      x_L=[x_L;zeros(n-length(x_L),1)]; 
   end
end
if length(x_U) < n
   if length(x_U)==1 
      x_U=x_U*ones(n,1);
   else
      x_U=[x_U;Inf*ones(n-length(x_U),1)]; 
   end
end

Prob.x_L = x_L;
Prob.x_U = x_U;
Prob.N   = n;

iFix=find(x_L==x_U);
iFree=find(isinf(x_L) & isinf(x_U));

if EMPTY
   Prob.x_0=zeros(n,1);
else
   Prob.x_0 = x;
end

% Safe-guard starting point
Prob.x_0    = max(Prob.x_L,min(Prob.x_0,Prob.x_U));

% ------------------------------------------------------------------------
function [x, B, Q, R, E, pRank, Stop] = GetXB(x, B, bl, bu, iFree, bFree,...
             bEqual, iFix, A, b_L, b_U, m, n0, EMPTY, rTol, epsRank, Prob)
% ------------------------------------------------------------------------

Stop=0;
[m,n]=size(A);

if EMPTY & isempty(B)
   x=bl;
   B=ones(n,1);
   m0=length(iFree);
   if m0 > 0
      x(iFree)=0;
      if m0 > m
         m0=m;
         B(iFree(1:m))=0;
      else
         B(iFree)=0;
      end
   end
   if sum(bFree) > 0
      x(n0+find(bFree))=0;
      B(n0+find(bFree))=0;
   end
   if ~isempty(iFix)
      x(iFix)=bl(iFix);
      B(iFix)=2;
   end
   B(isinf(x))=-1;
   x(isinf(x))=bu(isinf(x));

   if m0 > m
      % Too many degrees of freedom. Choose least norm solution
      fprintf('lpSimplex: Too many degrees of freedom.')
      fprintf(' %d free variables. m = %d\n',m0,m)
      % Compute the QR-decomposition
      [Q, R, E, pRank] = ComputeQR(A(:,1:n0), epsRank);
      x = tomsol(6, Q, R, E, pRank, b_L - A(:,1:n0)*x(1:n0));
      Stop=1;
      return;
   elseif length(iFree) == m
      % The base is OK
   else
      if m0==0
         % Just select the m slack variables as basis
         B(n0+1:n0+m)=0;
         m0=m;
      end
   end
elseif EMPTY
   % B is given
   if length(B) < n, B=[B;ones(n-length(B),1)]; end 
   if ~isempty(iFix)
      B(iFix)=2;
   end
   x=bl;
   x(find(B==-1))=bu(find(B==-1));
   if ~isempty(iFree)
      x(iFree)=0;
      B(iFree)=0;
   end
   m0=sum(B==0);
   if m0 > m
      fprintf(' Basis %d > m = %d\n',m0,m)
      Q=[]; R=[]; E=[]; pRank=0;
      Stop=2;
      return
   end
else
   % Use given initial x and possibly B
   if length(x) == n0
      x=[x;b_L];
      x(isinf(x))=bu(isinf(x));
      x(isinf(x))=0; % Free linear equations. Set slack to 0
   end
   if isempty(B) 
      B = ~((abs(x-bl) >  rTol .*max(1,abs(x))) &...
            (abs(x-bu) >  rTol .*max(1,abs(x))));
      
      B = double(B);

      B(find(abs(x-bu) <= rTol .*max(1,abs(x))))=-1;
      if ~isempty(iFix)
         x(iFix)=bl(iFix);
         B(iFix)=2;
      end
      if ~isempty(iFree)
         x(iFree)=0;
         B(iFree)=0;
      end

      m0=sum(B==0);
      %if m0 > m & size(A,1)~=0
      %   Q=[]; R=[]; E=[]; pRank=0;
      %   Stop=3;
      %   return
      %end

   else
      if length(B) < n 
         B=[B;ones(n-length(B),1)];
         B(n0+find(abs(x(n0+1:n0+m)-b_U) <= rTol(n0+1)))=-1;
      end
      if ~isempty(iFix)
         x(iFix)=bl(iFix);
         B(iFix)=2;
      end
      if ~isempty(iFree)
         x(iFree)=0;
         B(iFree)=0;
      end
      m0=sum(B==0);
      if m0 > m
         Q=[]; R=[]; E=[]; pRank=0;
         Stop=2;
         return
      end
   end
end
if m0 < m
   mq=sum(~bEqual & ~bFree);

   if m0+mq > m
      % First set the inequality slack variables as basis
      ix=find(~bEqual & ~bFree);
      B(n0+ix(1:m-m0))=0;
      m0=m;
   else
      B(n0+find(~bEqual & ~bFree))=0;
      m0=m0+mq;
   end
   if m0 < m
      ix=find(bEqual & ~bFree);
      B(n0+ix(1:m-m0))=0;
   end
end
   
% Compute the QR-decomposition
[Q, R, E, pRank] = ComputeQR(A(:,find(B==0)), epsRank);

% ------------------------------------------------------------------------
function [Alg,SolverAlgorithm] = SetAlg(Alg,SolverAlgorithm)
% ------------------------------------------------------------------------

if isempty(Alg), Alg=0; end

if Alg==1
   SolverAlgorithm=[SolverAlgorithm ' Bland''s rule.'];
elseif Alg==2
   SolverAlgorithm=[SolverAlgorithm ' Minimum reduced cost. Dantzig''s rule.'];
else
   Alg=0;
   SolverAlgorithm=[SolverAlgorithm ' Minimum reduced cost.'];
end
% ------------------------------------------------------------------------
function [x,B,Q,R,E,maxR,pRank] = HotRead(Prob,epsRank,PriLev)
% ------------------------------------------------------------------------


   if PriLev > 2, disp('Use Warm Start Basis'); end
   if ~isempty(Prob.QP.HotFile)
      load(Prob.QP.HotFile,'x','B','Q','R','E');
      if PriLev > 2
         fprintf('lpSimplex: Load Warm Start Basis File %s',Prob.QP.HotFile)
         fprintf('\n')
      end
   else
      if PriLev > 2
         fprintf('lpSimplex: Load Warm Start Basis From Prob structure')
         fprintf('\n')
      end
      x = Prob.QP.Hot.x;
      B = Prob.QP.Hot.B;
      Q = Prob.QP.Hot.Q;
      R = Prob.QP.Hot.R;
      E = Prob.QP.Hot.E;
   end
   maxR=max(abs(diag(R)));
   pRank = nnz(abs(diag(R)) >= maxR*epsRank);

% ------------------------------------------------------------------------
function [Ascale, A, D1, D2, b_L, b_U, x, x_L, x_U, c] = ScaleA(A, Ascale, ...
          iFix, iFree, b_L, b_U, bTol, x, x_L, x_U, c)
% ------------------------------------------------------------------------

[mA n]=size(A);
m=mA;

if isempty(Ascale), Ascale = 1; end

if mA == 0, Ascale = 0; end   % Safe guarding

if Ascale
   % Scaling of A
   iv = ones(n,1);
   if ~isempty(iFix)
      iv(iFix) = 0;
   end
   if ~isempty(iFree)
      iv(iFree) = 0;
   end
   ix = find(iv);

   if length(ix) > 1
      D1 = max(abs(A(:,ix))')';
      ix1 = find(D1==0);
      if ~isempty(D1),D1(ix1)=1; end

      ix1 = ~isinf(b_U);
      ix2 = ~isinf(b_L);
      if ~isempty(iFix)
         z = A(:,iFix)*x_U(iFix);
         D1(ix1) = max(D1(ix1),abs(b_U(ix1)-z(ix1)));
         D1(ix2) = max(D1(ix2),abs(b_L(ix2)-z(ix2)));
      else
         D1(ix1) = max(D1(ix1),abs(b_U(ix1)));
         D1(ix2) = max(D1(ix2),abs(b_L(ix2)));
      end

      for i = 1:m
          xx1 = (abs(A(i,:)) > bTol);
          minA = min(abs(A(i,iv & (abs(A(i,:)) > bTol)')));
          if isempty(minA)
             minA=1;
          end
          if ix1(i) & abs(b_U(i)) > bTol
             minA = min(minA,abs(b_U(i)));
          end
          if ix2(i) & abs(b_L(i)) > bTol
             minA = min(minA,abs(b_L(i)));
          end
          D1(i) = 1/sqrt(D1(i)*minA);
      end

%norm(A)
%plot(D1)
%pause

      A = diag(D1)*A;
      b_U(ix1) = D1(ix1).*b_U(ix1);
      b_L(ix2) = D1(ix2).*b_L(ix2);

      D2 = ones(n,1);

      D2(ix) = max(abs(A(:,ix)))';

      D2(D2==0)=1;

      for j = 1:length(ix)
          i = ix(j);
          minA = min(abs(A(abs(A(:,i)) > bTol,i)));
          if isempty(minA)
             minA=1;
          end
    
          D2(i) = 1/sqrt(D2(i)*minA);
      end
      A = A*diag(D2);
%norm(A)
%plot(D2)
%pause
      % Adjust bounds for the new variables
      x_L(ix) = x_L(ix) ./ D2(ix);
      x_U(ix) = x_U(ix) ./ D2(ix);
      if ~isempty(c)
         c(ix) = c(ix).*D2(ix);
      end
      if ~isempty(x)
         x(ix) = x(ix)./D2(ix);
      end
   else
      D1 = ones(m,1);
      D2 = ones(n,1);
   end
   D2=[D2;ones(m,1)];
end

%xprint(D1,'D1:')
%xprint(D2,'D2:')
%m
%rank(A)
%plot(svd(A))
%PriLev=1

% ------------------------------------------------------------------------
function [pOld,f_kOld,rSumOld,Cycle,doQR,qrFull,sinceQR,ixOld,ixL,ixU,y, ...
  Pen,pRank0,rAx,N_L,N_U,Iter] = InitVars(c,bl,bu,rTol,m,n0,B,pRank,PriLev)
% ------------------------------------------------------------------------

pOld   = 0;
f_kOld = NaN;
rSumOld= NaN;
Cycle  = 0;
doQR   = 0;
qrFull = 2000;
sinceQR= 0;

ixOld  = [];

n= length(c);

ixL    = zeros(n,1);
ixU    = zeros(n,1);
y      = zeros(m,1);

Pen = max(1000,2*max(abs(c))+1);
pRank0 = pRank;

if PriLev > 1
   fprintf('Rank at start %d. m = %d\n',pRank0,m);
end
rAx = ones(n-n0,1);
ix= find(~isinf(bl(n0+1:n)));
rAx(ix) = max(rAx(ix),bl(n0+ix));
ix= find(~isinf(bu(n0+1:n)));
rAx(ix) = max(rAx(ix),bu(n0+ix));

rAx=rTol(n0+1:n) .* rAx;

N_L = find(B==1);
N_U = find(B==-1);
Iter   = 0;

% ------------------------------------------------------------------------
function P1VarBug(P1Var,bl,bu,x,B,iB,Iter,n0,m,pRank,PriLev)
% ------------------------------------------------------------------------
if ~isempty(P1Var)
   disp('P1Var - ERROR in lpSimplex')
   [P1Var bl(P1Var) x(P1Var) bu(P1Var)]
   xprinti(iB,'iB:');
   xprinte(x(iB),'x(iB)');
   %pause
   if PriLev > 1
      fprintf('Iter %d. n %d. m %d.In iN %d ',Iter,n0,m,length(P1Var));
      fprintf('On Low %d. On Upp %d pRank %d\n',... 
               length(find(B==0 & x==bl)),  ...
               length(find(B==0 & x==bu)), pRank  ...
      );
   end
end

% ------------------------------------------------------------------------
function ResultLP = SaveResult(flag, Ascale, x, y, D1, D2, B, Q, R, E,...
   Iter, n0, m, f_k, cHat, N_U, pRank, x_L, x_U, ResultLP)
% ------------------------------------------------------------------------

ResultLP.x_k=xScale(x, Ascale, D2, n0);
if Ascale
   ResultLP.y_k=y./D1;
else
   ResultLP.y_k=y;
end
ResultLP.FuncEv   = Iter;
ResultLP.ConstrEv = Iter;
if flag ~= 2000
   ResultLP.QP.B=B(1:min(length(B),n0));
   ResultLP.QP.Q=Q;
   ResultLP.QP.R=R;
   ResultLP.QP.E=E;
   ResultLP.QP.pRank=pRank;
end
if ~isempty(N_U) & ~isempty(cHat)
   % Lagrange multipliers for bounds is reduced costs
   cHat(N_U)=-cHat(N_U);
end
if flag==4
   ResultLP.v_k=zeros(n0,1);
elseif flag==3
   ResultLP.v_k=[];
elseif flag==1000
   ResultLP.v_k=zeros(n0,1);
   flag=0;
   ResultLP.QP.BFull=B;
elseif flag==2000
   ResultLP.v_k=zeros(n0,1);
   flag=0;
else
   ResultLP.v_k=[cHat(1:n0);y];
   ResultLP.QP.BFull=B;
end
ResultLP.ExitFlag=flag;
ResultLP.ExitText=ExitText(flag);
ResultLP.f_k=f_k;
%ResultLP.xState=(x(1:n0)==x_L)+2*(x(1:n0)==x_U);
ResultLP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );

% ------------------------------------------------------------------------
function x = xUpdate(x, q, bl, bu, alpha, Implicit, qVar, LOW)
% ------------------------------------------------------------------------
if alpha > 0 & ~Implicit
   if LOW
      x(q) = bl(q) + alpha;
   else
      x(q) = bu(q) - alpha;
   end
elseif alpha > 0 & qVar~=0
   if LOW
      x(q) = x(q) + alpha;
   else
      x(q) = x(q) - alpha;
   end
end
% ------------------------------------------------------------------------
function [iB,N_L,N_U,iN] = Bsets(B)
% ------------------------------------------------------------------------

iB  = find(B==0);
if nargout < 2, return, end
N_L = find(B==1);
N_U = find(B==-1);
if nargout < 4, return, end
iN  = find(B==1 | B==-1);

% ------------------------------------------------------------------------
function [x] = xScale(x, Ascale, D2, n0)
% ------------------------------------------------------------------------
if Ascale
   x=x(1:n0).*D2(1:n0);
else
   x=x(1:n0);
end

% ------------------------------
function Text = ExitText(Inform)
% ------------------------------

switch  Inform
   case 0
     Text = 'Optimal solution found';
   case 1
     Text = 'Maximal number of iterations reached';
   case 2
     Text = 'Unbounded feasible region';
   case 4
     Text = 'Infeasible problem wrt bounds or constraints';
   case 5
     Text = 'Too many active variables in given initial point';
   case 6
     Text = 'No feasible point found with Phase 1';
   case 10
     Text = 'Errors in input parameters';
   case 11
     Text = 'Illegal initial x as input';
   otherwise
     Text = [];
end

% MODIFICATION LOG:
%
% 990809  hkh   Revised for 2.0
% 990821  hkh   Completely rewritten, slacks for constraints, new Phase 1.
% 990901  hkh   FullRank now works. Tested on knapsack problems.
% 990913  hkh   Safeguard against inf and nan in x_0
% 000916  hkh   Add text about convergence to ExitText
% 020820  ango  Change for Matlab 6.5 logical handling
% 020821  hkh   Apostrophes missing in save statement
% 030416  ango  Further 6.5 logical handling fixes (FullRank)
% 030729  ango  Logical handling again
% 040111  hkh   Add call to inisolve and endSolve
% 040923  hkh   format compact long changed to format compact, not in R13-14
% 041222  med   Safeguard added for x_0
% 050117  med   mlint revision
% 051215  hkh   Exit text missing for case 4, infeasible problem
% 051215  hkh   Inform value for ExitFlag=4 not documented
% 060816  med   Removed printing by default
% 080606  med   Switched to iniSolveMini