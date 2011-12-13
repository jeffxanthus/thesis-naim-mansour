% Phase II - Dual simplex algorithm.
%
% Given LP on the standard form, but with general lower and upper bounds:
%       min c' * x , A x = b_U, x_L <= x <= x_U
%
%
% Solve dual LP:   max  b'*y, subject to  A'*y <= c, y urs
%
% Note! A known dual feasible initial basis is needed in Prob.QP.B
%
% Upper bounds (not Inf)  are added as extra linear constraints in A.
% Lower bounds (not zero) are added as extra linear constraints in A.
% Fixed variables are treated explicitely.
%
% Note that the TOMLAB routine cpTransf.m can be used to rewrite
% an LP problem on another form to this standard form
%
% Implementation of algorithm Dual simplex method,
% Handbooks in Operations Research and Management Science, Chapter II:
% Goldfarb and Todd: Linear Programming, pages 105-106
%
% Generalized with  QR factorization, and numerical safeguarding
%
% function Result = DualSolve(Prob);
%
% INPUT:
% In the Prob structure use:
%
% x_L  Lower bound on x.
%      Assumed to be 0, if not Prob.x_L(i) == Prob.x_U(i) (fixed vars)
% x_U  Upper bound on x.
%
% b_L  Lower bound on linear constraints. Not used.
% b_U  Upper bounds, Assume equality constraints Ax == b_U
%
% A,b_U,c   A in R^(m x n). b_U in R^m. x in R^n.
%
% x_0  Starting value x_0
%
% QP.B Active set B at start.
%           1  = Include variable x(i) is in basic set.
%           0  = Variable x(i) is set on its lower bound 0
%           -1 = Variable x(i) is set on its upper bound (NOT ALLOWED)
%
% QP.y Dual variables y (Lagrange multipliers for linear constraints)
%
% Starting value set {x_0,B,y} is given in the following three ways:
%
% {x_0,B,y} All defined
%
% {x_0,B}   Estimate y from A(:,iB)'y=c(iB), where basis iB=find(B==1)
%
% {B}       x_0 is estimated from A(:,iB) x = b_U. y from A(:,iB)'y=c(iB)
%
% PriLevOpt Printing level:
%           =0 No output; >0 Convergence results;
%           >1 Output every iteration  >2 Output each step in the simplex alg
%
% Prob.Solver.Alg  Rule to select new variables:
%           = 0 Minimum Reduced Cost (DEFAULT), sort variables increasing.
%           = 1 Bland's Anti-cycling Rule
%           = 2 Minimum Reduced Cost, Dantzig's rule
%
% Fields used in Prob.optParam
%
% MaxIter   Maximal number of iterations. max(10*dim(x),100) is DEFAULT.
% eps_f     Tolerance used to test convergence on the reduced costs
% eps_Rank  Rank tolerance
% xTol      Tolerance to judge if x-values are close
% bTol      Tolerance for linear constraints
% wait      Wait flag, pause each iteration if set true
%
% OUTPUT:
% Result structure with the following fields defined:
% Field Internal
% name  name
% x_k    x    Primal solution x
% y_k    y    Dual solution y. Lagrangian multipliers for the linear constraints
% QP.B   B    Optimal set. See QP.B as INPUT.
% QP.DualLimit Limit value for b'*y. If b'y >= QP.DualLimit, Stop.
% ExitFlag    Exit flag
%             == 0  => OK
%             == 1  => Maximal number of iterations reached. No bfs found.
%             == 2  => Infeasible dual problem
%             == 4  => Too many active variables in initial point given.
%             == 6  => No dual feasible starting point found.
%             == 7  => Illegal step length due to numerical difficulties.
%                      Should not occur.
%             == 8  => Convergence when f_k >= QP.DualLimit. No further
%                      progress is wanted. Used by MIP codes.
%             == 9  => x_L(i) > x_U(i) + xTol for some i. No solution exists
% ExitText    Text string giving ExitFlag and Inform information
%
%             Some reduced cost is negative.
% f_k         f_hat=b'y at optimum y (or last iterate y if no convergence)
% A           Matrix A in standard LP formulation.
% b           Vector b in standard LP formulation.
% c           Vector c in standard LP formulation.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1994-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.1.0$
% Written Nov 10, 1994.   Last modified Dec 16, 2005.

function ResultDLP = DualSolve(Prob)

if nargin < 1 
   error('DualSolve needs one input parameter Prob');
end

solvType=checkType('lp');

Prob=ProbCheck(Prob,'DualSolve',solvType,8);

ResultDLP                 = ResultDef(Prob);
ResultDLP.Solver          = 'DualSolve';
ResultDLP.Prob            = Prob;
ResultDLP.SolverAlgorithm = 'Dual Simplex method.';

x   = Prob.x_0(:);  
x_L = Prob.x_L(:);
x_U = Prob.x_U(:);

n = length(x);                  % # of variables

if n == 0
   n = max(length(x_L),length(x_U));
end

if isempty(x), x = zeros(n,1); end

if isempty(x_U),x_U= Inf*ones(n,1); end
if isempty(x_L),x_L= zeros(n,1); end

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
% Safe-guard starting point
x = max(x_L(1:n),min(x,x_U(1:n)));

Prob.x_0 = x;
Prob.x_L = x_L;
Prob.x_U = x_U;
Prob.N   = n;

x_L=max(0,x_L);

xEqual=eq(x_L,x_U);

optParam=Prob.optParam;
xTol     = optParam.xTol;
bTol     = optParam.bTol;

if any(x_L > x_U + xTol)
   ResultDLP.f_k=NaN;
   ResultDLP.QP.B=[];
   ResultDLP.Iter=0;
   ResultDLP.ExitFlag=9;
   ResultDLP.ExitText=ExitText(9);
   return
end

A   = Prob.A;
b_U = Prob.b_U(:);
b_L = b_U;
Prob.b_L = b_U;
c   = Prob.QP.c(:);

[m, n] = size(A);

Ascale   = Prob.QP.Ascale;

if isempty(Ascale), Ascale = 1; end

if Ascale
   % Scaling of A

   D1 = max(abs(A)')';

   ix1 = ~isinf(b_U);
   ix2 = ~isinf(b_L);

   D1(ix1) = max(D1(ix1),abs(b_U(ix1)));
   D1(ix2) = max(D1(ix2),abs(b_L(ix2)));

   for i = 1:m
       minA = min(abs(A(i,abs(A(i,:)) > bTol)));
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

   %norm(A)

   D2 = max(abs(A))';
   D2(D2==0)=1;

   for i = 1:n
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
   x_L = x_L ./ D2;
   x_U = x_U ./ D2;
   if ~isempty(c)
      c = c.*D2;
   end
   if ~isempty(x)
      x = x./D2;
   end
   D2=[D2;ones(m,1)];
end


if isempty(c)
    c=zeros(n,1);
end

B = Prob.QP.B(:);
if length(B)~=n 
   error('Prob.QP.B must contain a correct basis');
end
if all(x==0)
   xprinte(B,'B:');
end
y = Prob.QP.y(:);

DualLimit = Prob.QP.DualLimit;
if isempty(DualLimit), DualLimit=Inf; end

if isempty(Prob.P), Prob.P=1; end

%find(B==-1 & abs(x) < xTol & abs(x_U) < xTol)
%find(B==-1 & abs(x) < xTol )
%B(find(B==-1 & abs(x) < xTol & abs(x_U) < xTol))=0;
%B(find(B==-1 & abs(x) < xTol))=0;
%B(find(B==-1 & abs(x) > xTol & abs(x-x_U) > xTol))=1;

% Fixed variables

iEQ=find(xEqual & abs(x - x_U) < xTol);
iEE=find(xEqual & abs(x - x_U) >= xTol);
mEQ=length(iEQ);
mEE=length(iEE);
mEQ0=0;
if mEQ > 0
   mEQ0=sum(x_U(iEQ)==0);
   x(iEQ) = x_U(iEQ);
   b_UEQ = b_U - A(:,iEQ)*x(iEQ);
else
   b_UEQ = b_U;
   yEQ=[];
end

iU=find((~xEqual) & (~isinf(x_U)));
iL=find((~xEqual) & (x_L~=0) & (~isinf(x_L)));

m0=m;
n0=n;
mU=length(iU);
mL=length(iL);

% Add upper bounds (not equal to Inf) as linear constraints
if mU > 0
   b_U = [b_U;x_U(iU)];
   b_UEQ = [b_UEQ;x_U(iU)];
   A = [A zeros(m,mU); full(sparse(1:mU,iU,ones(mU,1),mU,n)) eye(mU,mU)];
   n = n + mU;
   m = m + mU;
   x = [x;x_U(iU)-x(iU)];
   c = [c;zeros(mU,1)];
   y = [];
end
% Add lower bounds (not equal to zero) as linear constraints
if mL > 0
   b_U = [b_U;x_L(iL)];
   b_UEQ = [b_UEQ;x_L(iL)];
   A = [A zeros(m,mL); full(sparse(1:mL,iL,ones(mL,1),mL,n)) -eye(mL,mL)];
   n = n + mL;
   m = m + mL;
   x = [x;x(iL)-x_L(iL)];
   c = [c;zeros(mL,1)];
   y = [];
end
% Add fixed variables, where x is not on the fixed value
if mEE > 0
   b_U = [b_U;x_L(iEE)];
   b_UEQ = [b_UEQ;x_L(iEE)];
   iz = find(x(iEE) < x_U(iEE));
   % To guarantee dual feasibility, slack s must always be negative.
   z  = ones(mEE,1);  % x + s = x_U, if x_U < x
   z(find(x(iEE) < x_U(iEE)))=-1;   % x - s = x_U, if x_U > x
   x = [x;-abs(x(iEE)-x_L(iEE))];
   A = [A zeros(m,mEE);full(sparse(1:mEE,iEE,ones(mEE,1),mEE,n)) ...
                       full(sparse(1:mEE,1:mEE,z,mEE,mEE))];
   n = n + mEE;
   m = m + mEE;
   c = [c;zeros(mEE,1)];
   y = [];
end

wait     = optParam.wait;     % Pause after printout if true
PriLev   = Prob.PriLevOpt;       % Print level
epsRank  = optParam.eps_Rank; % Rank test tolerance
MaxIter  = optParam.MaxIter;  % Maximal number of iterations
eps_f    = optParam.eps_f;    % Reduced cost convergence tolerance

if PriLev > 3
   xprinti(b_U(m0+mU+1:m),'b')
   mPrint(A(m0+mU+1:m,:))
   format long
   x(n0+1:n)
   mU
   if wait, pause; end
end

if isempty(xTol), xTol = 100*eps; end

DEBUG=0;


Alg=Prob.Solver.Alg;

if isempty(Alg), Alg=0; end

Prob.Solver.Alg=Alg;

if Alg==0
   ResultDLP.SolverAlgorithm=...
          [ResultDLP.SolverAlgorithm ' Minimum reduced cost.'];
elseif Alg==1
   ResultDLP.SolverAlgorithm=...
          [ResultDLP.SolverAlgorithm ' Bland''s rule.'];
elseif Alg==2
   ResultDLP.SolverAlgorithm=...
          [ResultDLP.SolverAlgorithm ' Minimum reduced cost. Dantzig''s rule.'];
end

if PriLev > 2
   fprintf('Run %s',ResultDLP.SolverAlgorithm);
   fprintf('\n');
end

IndexRule = Prob.Solver.Alg;

B=ones(n0,1);
B(abs(x(1:n0)) > xTol)=0;

if mU > 0
   B = [B;ones(mU,1)];
   ixx=find(B(iU)==1 | abs(x(n0+[1:mU])) >= xTol);
   B(n0+ixx)=0;
   %B = [B;~(abs(x(n0+[1:mU])) >= xTol)];
end
if mL > 0
   B = [B;~(abs(x(n0+mU+[1:mL])) >= xTol)];
end
if mEE > 0
   B = [B;zeros(mEE,1)];
   B(iEE)=0;
end

if mEQ > 0
   B(iEQ) = 2;
end

if any(b_L < -xTol)
   if PriLev > 0
      disp('Negative b_L')
      xprinte(b_L,'b_L');
   end
end

mm=sum(B==0);
if mm < m
   % Try to increase
   ix=find(B(1:n0)==1 & (Prob.QP.B==0 | Prob.QP.B==-1));
   if length(ix)==m-mm
      B(ix)=0;
   else
      [iy iz]=sort(abs(x(ix)));
      nn=min(length(ix),m-mm);
      B(ix(iz(1:nn)))=0;
   end
end


iN=find(B==1);
iB=find(B==0);

if length(iB) < m 
   if PriLev > 0
      fprintf('Problem with initial B. Too few variables. ');
      fprintf('Problem %d\n',Prob.P);
      fprintf('Basis is %d. Should be m = %d\n',length(iB),m);
   end
elseif length(iB) > m

   if PriLev > 0
      disp('Problem with initial B. Too many variables');
      fprintf('Problem %d\n',Prob.P);
      fprintf('Basis is %d. Should be m = %d\n',length(iB),m);
   end

if 0
   ixx=find(B(n0+1:n0+mU+mL)==0 & abs(x(n0+1:n0+mU+mL)) < xTol);

   % Delete slack variables from the base to make it have size m
   if PriLev > 0
      fprintf('Delete slack variables from the base, making it have size m')
      fprintf(' %d from %d vars\n',m,length(iB));
   end

   if ~isempty(ixx)
      k=length(iB)-m;
      B(n0+ixx(1:k))=1;
   else
      ixx=find(B==1 & abs(x) < xTol);
      if ~isempty(ixx)
         k=min(length(ixx),length(iB)-m);
         B(n0+ixx(1:k))=1;
      end
   end
   iN=find(B==1);
   iB=find(B==0);
end
end

rx=x(1:n0)-Prob.x_0(1:n0);

% Make QR factorization

[Q, R, E, pRank] = ComputeQR(A(:,iB), epsRank);

if length(iB) > m
   % Remove variables from basis, based on the result from the QR
   [P,jE]=find(E);
   % Remove extra unnecessary columns
   for i=length(iB):-1:m+1
       % Drop number i. Renumber P
       j = P(i); 
       k = iB(j);
       x(k)=0;
       B(k)=1;
       P(P > j)=P(P > j) - 1;
       iB=find(B==0);
   end
   P=P(1:m);
   R=R(:,1:m);
   iN=find(B==1);
   E=sparse(P(:),[1:m]',ones(m,1),m,m);
end

if isempty(iN)
   x(iB) = tomsol(6, Q, R, E, pRank, b_UEQ);
else
   x(iB) = tomsol(6, Q, R, E, pRank, b_UEQ-A(:,iN)*x(iN));
end

if sum(abs(B)==0) < m | pRank < m
   pRank0=pRank;
   % Generate a full rank basis, if possible
   [B, Q, R, E, pRank, DelVar, NewVar] = FullRank(B, A, Q, R, E, pRank,...
       m, epsRank);
   if PriLev > 1 & pRank0 < pRank
      fprintf('FullRank: Increased rank from %d to %d.',pRank0,pRank);
      fprintf(' Rows m = %d.\n',m);
   end
   x(DelVar)=0;
   x(NewVar)=0;

   %  r_k=b_U-A*x;
   %  % Try a least squares solution
   %  x(NewVar)=A(:,NewVar)\r_k;
   %if norm(b_U-A*x) < 1E-12 * max(1,norm(b_U))

   if pRank==m
%NewVar
      for i=1:length(NewVar)
          k=NewVar(i);
          r_k=b_U-A*x;
          ix = find(A(:,k) ~= 0);
          New = r_k(ix)./A(ix,k);
%xprint(New,'New');
          [maxNew j]=max(New);
          if maxNew==0
             x(k)=max(r_k)/A(ix(j),k);
%x(k)
          else
             x(k)=New(j);
          end
%pause
      end
   else
      % We give up and call lpSimplex instead
      %disp('**** NOT FULL RANK, call lpSimplex **************')
      Prob.x_0=[];
      Prob.B  =[];
      ResultDLP=lpSimplex(Prob);
      return
   end
      

   iB  = find(B==0);
   iN = find(B==1);
   %N_U = find(B==-1);
end

if ~isempty(y)
   r=c-A'*y;
   rSum=sum(abs(r));
   if rSum > 1E-12 * max(1,sum(abs(c)))
      if PriLev > 0
         disp('DualSolve: Must recompute dual parameters y!!! ');
         fprintf('Given y is not feasible. sum(|c-A''b|) = %20.10e\n',rSum);
      end
      y=[];
   end
end

if isempty(y)
   % Step 0
   % ----- Solve the linear equation  A_B^T * y = c_B for the dual vector y
   % y = A(:,iB)' \ c(iB);
   % QR - factorization already computed: [Q R E] = qr(full(A(:,iB)));
   y = tomsol(7,Q, R, E, pRank, c(iB));

   if PriLev > 4
      fprintf('Shadow prices:\n');
      xprinte(y,'y: ');
   end
end
if mEQ > 0
   yEQ=c(iEQ)-A(:,iEQ)'*y;
end


if PriLev > 2
   xprinti(iB,'iB:');
   xprinti(iN,'iN:');
end

% ----- Compute reduced costs 

cHat = c - A' * y;

if PriLev > 4
   fprintf('Reduced costs:\n');
   xprinte(cHat(iN),'cHat:');
end
if any(abs(cHat(iB)) > 1E-5)
   if PriLev > 4
      fprintf('Reduced costs for basis should be 0. ');
      fprintf('Max is %e\n', max(abs(cHat(iB))));
      disp('Recompute QR');
   end
   doQR=1;
end

if any(cHat(iN) < 0)
   if any(cHat(iN) < -eps_f)
      if PriLev > 0
         disp('Error in Dual Simplex. Shadow prices y NOT DUAL FEASIBLE!!!')
         fprintf('Shadow prices:\n');
         xprinte(y,'y:    ');
         fprintf('Reduced costs:\n');
         xprinte(cHat(iN),'cHat:');
      end
      ResultDLP.x_k=xScale(x, Ascale, D2, n0);
      ResultDLP.y_k=[];
      ResultDLP.f_k=NaN;
      B(iEQ)=0;
      ResultDLP.QP.B=B(1:n0);
      %ResultDLP.QP.B=[];
      ResultDLP.Iter=0;
      ResultDLP.ExitFlag=6;
      ResultDLP.ExitText=ExitText(6);

      % EMERGENCY! Call lpSimplex
      if PriLev > 0
         disp('*** NOT DUAL FEASIBLE! Call lpSimplex instead ***')
      end
      Prob.x_0=[];
      Prob.B  =[];
      ResultDLP=lpSimplex(Prob);
      return
   else  % Avoid problems due to numerical rounding
      %disp('AVOID NUMERICAL PROBLEMS IN DualSolve')
      j=find(cHat(iN) > -eps_f & cHat(iN) < 0);
      cHat(iN(j))=zeros(length(j),1);
   end
end

ResultDLP.y_k=y(1:m0);
sinceQR=0;
doQR=0;

pRank0 = pRank;

if PriLev > 1
   fprintf('Rank at start %d. m = %d\n',pRank0,m);
end
% NOTE! Change m if A is rank deficient.
m=pRank0;

bNorm=norm(b_U);

for Iter = 1:MaxIter
    f_k=b_U'*y;
    r=A*x-b_U;

    if any(abs(r) > 1E-10*max(1,abs(b_U)))  % Try to restore feasibility
       if PriLev > 0
          disp('Problem with feasibility')
          rSum=sum(abs(r));
          rSum
       end
    end
    if mEQ > 0
       f_k=f_k+x_U(iEQ)'*yEQ;
    end
    if PriLev > 1
       fprintf('\n--------------------------------------------------\n')
       fprintf('DUAL Simplex, Problem %d. Iteration: %5.0f\n',Prob.P,Iter)
       fprintf('--------------------------------------------------\n')
       fprintf('Max LP b''*y =     %30.16f\n',f_k);
       % cTx always equal to b'*y
       cTx=c'*x;
       fprintf('Min LP c''*x =     %30.16f\n',cTx);
       fprintf('Duality gap =     %30.16f\n',f_k-cTx);
       fprintf('n0 %d m0 %d n %d m %d pRank %d ',n0,m0,n,m,pRank);
       fprintf('iB %d mEQ %d mU %d ',length(iB),mEQ,mU);
       fprintf('mEQ0 %d mEE %d mL %d\n',mEQ0,mEE,mL);
       fprintf('\n');

       if PriLev > 2
          fprintf('Base B:\n');
          xprinti(iB,'iB:');
          fprintf('Non-basic set N:\n');
          xprinti(iN,'iN:');
          xprint(x,'x:    ')
       end
       if PriLev > 3
          xprint(cHat,'cHat:');
          xprint(y,'y:    ');
       end

       if wait
          disp('Press any key to continue ...')
          pause 
       end
    end
    ResultDLP.x_k=xScale(x, Ascale, D2, n0);
    ResultDLP.y_k=yScale(y, Ascale, D1, m0);
    ResultDLP.f_k=f_k;
    ResultDLP.QP.B=B(1:n0);
    ResultDLP.Iter=Iter;
  
    % ===== Step 1: Convergence test. Check primal feasibility =====
    if all(x >= -xTol)
       if PriLev > 0
          fprintf('DualSolve: Problem %d. x is primal feasible.',Prob.P)
          fprintf(' Convergence!\n')
       end
       if PriLev == 2
          xprint(x,'x:    ')
       end
       if any(x) < 0 % Fix small x to exactly 0
         if PriLev > 0
          disp('x is primal feasible, but with very small negative values.')
          disp('Assume Convergence!!!')
          xprinte(x,'x:    ')
         end
          j=find(x < 0);
          x(j)=zeros(length(j),1);
          if PriLev > 0
             fprintf('Max LP b''*y =     %30.16f\n',f_k);
             % cTx always equal to b'*y
             cTx=c'*x;
             fprintf('Min LP c''*x =     %30.16f\n',cTx);
             fprintf('Duality gap =     %30.16f\n',f_k-cTx);
          end
       end
       if PriLev > 1
          xprint(iB,'Base B:         ',' %4.0f');
          xprint(iN,'Non-basic set N:',' %4.0f');
       end
       if PriLev > 3
          xprint(cHat,'cHat:');
          xprint(y,'y:    ');
       end
       % Compute full Lagrange multiplier vector
       v_k=cHat(1:n0);
       v_k(iU)=y(m0+1:m0+mU);
       v_k(iL)=y(m0+mU+1:m0+mU+mL);
       v_k(iEE)=y(m0+mU+mL+1:m0+mU+mL+mEE);
       v_k(iEQ)=yEQ;
       ResultDLP.x_k=xScale(x, Ascale, D2, n0);
       ResultDLP.y_k=yScale(y, Ascale, D1, m0);
       ResultDLP.v_k=[v_k;y(1:m0)];
       ResultDLP.f_k=f_k;
       B(iEQ)=0;
       if mU > 0
          iz0=find((B(n0+[1:mU])==1) & (abs(x(iU)-x_U(iU)) < xTol)  ...
                   & (x_U(iU)>xTol) & (~isinf(x_U(iU))));
%xprinti(iz0,'iz0')
%xprinti(B(iU(iz0)),'B-0');
%xprinte(x(iU(iz0)),'x-1');
%xprinte(x_U(iU(iz0)),'xU');
          B(iU(iz0))=-1;
       end
       ResultDLP.QP.B=B(1:n0);
       ResultDLP.QP.Q=Q;
       ResultDLP.QP.R=R;
       ResultDLP.QP.E=E;
       ResultDLP.QP.pRank=pRank;
       ResultDLP.Iter=Iter;
       ResultDLP.ExitFlag=0;
       ResultDLP.ExitText=ExitText(0);
       return;
    end
    if f_k >= DualLimit
       if PriLev > 0
          disp('DualSolve: DualLimit exceeded. Stop iterations!')
       end
       if PriLev > 1
          xprint(x,'x:    ')
       end
       ResultDLP.x_k=xScale(x, Ascale, D2, n0);
       ResultDLP.y_k=yScale(y, Ascale, D1, m0);
       ResultDLP.f_k=f_k;
       B(iEQ)=0;
       ResultDLP.QP.B=B(1:n0);
       ResultDLP.Iter=Iter;
       ResultDLP.ExitFlag=8;
       ResultDLP.ExitText=ExitText(8);
       return;
    end

    % ===== Step 2:  Choose index for variable x to exclude from basis =====
    if IndexRule == 1    		
       % Choose index using Blands rule
       pHat = find(x(iB) < -1E-10);
       if isempty(pHat), pHat = find(x(iB) < 0); end
       pHat = pHat(1);
    else                   		
       % Choose index using Minimal x value
       [min_x pHat] =  min(x(iB));
    end
    p = iB(pHat);
    % ===== Step 3:  Check infeasibility
    e_p = zeros(m,1);
    e_p(pHat) = 1;
    % ----- Solve the linear equation  A_B' * u = e_pHat. Row pHat in A_B^(-1)
    % QR - factorization already computed: [Q R E] = qr(full(A(:,iB)));
%if DEBUG
%    u1 = A(:, iB)' \ e_p;
%end
   %===== Step 3:  Compute shadow prices =========================
   % y = A(:,iB)' \ c(iB);
   % QR - factorization already computed: [Q R E] = qr(full(A(:,iB)));

   u = tomsol(7,Q, R, E, pRank, e_p);

%if DEBUG
%   xprinte(u,'u: ');
%   xprinte(u1,'u1:');
%   xprinte(u1-u,'uD:');
%   sum(abs(u1-u))
%end

    v = A(:, iN)' * u;
    if PriLev > 3
       disp('u = A_B(-1) * e_pHat --- search direction in dual variables y:')
       xprinte(u,'u:    ')
       disp('v = A_N'' * u          --- search direction in reduced costs:')
       xprinte(v,'v:    ')
    end
    if PriLev > 2
       fprintf('Excluded variable x_p:           ')
       fprintf('p = %5.0f, pHat = %5.0f \n',p,pHat)
    end
    %if all(v >= -eps_f)
    if all(v >= -100*eps)
       if PriLev > 0
          if PriLev == 1 | PriLev == 2 | PriLev == 3 % Otherwise already printed
             fprintf('v = A_N'' * u          ')
             fprintf('--- search direction in reduced costs:\n')
             xprinte(v,'v:    ')
          end
          disp('Infeasible DUAL problem!!!');
       end
       ResultDLP.x_k=xScale(x, Ascale, D2, n0);
       ResultDLP.y_k=yScale(y, Ascale, D1, m0);
       ResultDLP.f_k=f_k;
       B(iEQ)=0;
       ResultDLP.QP.B=B(1:n0);
       ResultDLP.Iter=Iter;
       ResultDLP.ExitFlag=2;
       ResultDLP.ExitText=ExitText(2);
       return;
    end
    % ===== Step 4:  Determine non basic variable x_q to enter in base
    gamma = Inf;
    for k = 1:length(iN)
        if v(k) < -1E-10
           gamma_k = -cHat(iN(k)) / v(k);
           if gamma_k < gamma
              gamma = gamma_k;
              qHat = k;
           end
        end
    end
    if isinf(gamma)
       % TRY TO AVOID NUMERICAL PROBLEMS DUE TO SMALL v
       for k = 1:length(iN)
           if v(k) < -1E-15
              gamma_k = -cHat(iN(k)) / v(k);
              if gamma_k < gamma
                 gamma = gamma_k;
                 qHat = k;
              end
           end
       end
    end
    if gamma == Inf
       % THIS SHOULD NOT HAPPEN!!!
       disp('GAMMA == INF')
       if PriLev > -1
          fprintf('DualSolve: ERROR! This should not occur! \n')
          fprintf('No v < 0 ???? \n')
          xprinte(x,'x:    ')
          xprinte(y,'y:    ')
          xprinte(u,'u:    ')
          fprintf('Quote between cHat and v too high.\n')
          xprinte(cHat(iN),'cHat:')
          xprinte(v,'v:    ')
          % cTx always equal to b'*y
          cTx=c'*x;
          fprintf('Min LP c''*x =     %30.16f\n',cTx);
          fprintf('Duality gap =     %30.16f\n',f_k-cTx);
       end
       ResultDLP.x_k=xScale(x, Ascale, D2, n0);
       ResultDLP.y_k=yScale(y, Ascale, D1, m0);
       ResultDLP.f_k=f_k;
       B(iEQ)=1;
       ResultDLP.QP.B=B(1:n0);
       ResultDLP.Iter=Iter;
       ResultDLP.ExitFlag=7;
       ResultDLP.ExitText=ExitText(7);
       return;
    end
    q = iN(qHat);

    if PriLev > 2
       fprintf('Included variable x_q:           ')
       fprintf('q = %5.0f, qHat = %5.0f \n',q,qHat)
       fprintf('Step in dual   parameters is gamma = %20.10f \n',gamma)
    end
    if gamma > 1E10 % Numerical difficulties.
       fprintf('Numerical difficulties computing gamma %20.10e\n',gamma)
       if PriLev > 2
          xprinte(x,'x:    ')
          xprinte(y,'y:    ')
          xprinte(u,'u:    ')
       end
       fprintf('Quote between cHat and v too high. Exit DualSolve.\n')
       if PriLev > 2
          xprinte(cHat(iN),'cHat:')
          xprinte(v,'v:    ')
       end
       if PriLev > 1
          % cTx always equal to b'*y
          cTx=c'*x;
          fprintf('Min LP c''*x =     %30.16f\n',cTx);
          fprintf('Duality gap =     %30.16f\n',f_k-cTx);
       end
       ResultDLP.x_k=xScale(x, Ascale, D2, n0);
       ResultDLP.y_k=yScale(y, Ascale, D1, m0);
       ResultDLP.f_k=f_k;
       B(iEQ)=0;
       ResultDLP.QP.B=B(1:n0);
       ResultDLP.Iter=Iter;
       ResultDLP.ExitFlag=4;
       ResultDLP.ExitText=ExitText(4);
       return
    end

    % ===== Step 5:  Update reduced costs and y
    cHat(iN) = cHat(iN) + gamma * v;
    cHat(p) = gamma;
    cHat(q) = 0;    % Set to identically 0 to avoid round off difficulties 
    y = y - gamma * u;
    if mEQ > 0
       yEQ=c(iEQ)-A(:,iEQ)'*y;
    end

    % ===== Step 6:  Update solution x and basis B
%if DEBUG
%    d1 = A(:,iB) \ (-A(:, q));
%end

   d = tomsol(6, Q, R, E, pRank, -A(:,q));

   dErr=sum(abs(A(:,iB)*d+A(:,q)));
   if PriLev > 1
      fprintf('|A(:,iB)*d +- A(:,q)| = %23.10e\n',dErr);
   end
   if PriLev > 4
      xprinte(d,'d:');
   end

   if dErr > 1E-8*bNorm
      if PriLev > 1 | dErr > 0.1*bNorm
         disp('Residual for search direction too large. Force QR rebuild')
         fprintf('Feasibility? 0 = |A*x - Is | = %23.10e\n',norm(A*x));
         %% We give up and call lpSimplex instead
         %disp('**** Numerical ill-conditioning, call lpSimplex **************')
         %Prob.x_0=[];
         %Prob.B  =[];
         %ResultDLP=lpSimplex(Prob);
         %return
      end
      doQR=1;
   end

%if DEBUG
%   xprinte(d,'d: ');
%   xprinte(d1,'d1:');
%   xprinte(d1-d,'dD:');
%   sum(abs(d1-d))
%   pause
%end

   alpha = x(p)/v(qHat);

   if PriLev > 2
      fprintf('Step in primal parameters is alpha = %20.10f \n',alpha)
   end
   if PriLev > 3
      fprintf('d --- Search direction in x:\n')
      xprint(d,'d:    ')
   end
   x(q) = alpha;
   x(iB) = x(iB) + alpha * d;
   x(p) = 0;        % Set to identically 0 to avoid round off difficulties 
   B(p) = 1;
   B(q) = 0;

   if IndexRule <=1
      iB = find(B==0);
      iN = find(B==1);

   else                     % Dantzig's rule
      iN(qHat) = p;
      iB(pHat) = q;
   end
   if doQR > 0 
      if doQR > 0
         if PriLev > 1
            fprintf('FORCED QR!!!\n');
         end
         doQR=0;
      else
         if PriLev > 1
            fprintf('Rebuild QR!!! Iter %d\n',Iter);
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

      [Q,R,P,pRank,maxR,absR] = InsertQR(Q, R, P, pRank, k, A(:,q), epsRank);


%if any(sort(P)~=[1:m]')
%   disp('DualSolve: Serious P error. Something is totally wrong.');
%   %disp('Bad input data, or some logical program error.');
%   xprinti(P,'P:');
%   pause
%end
      if pRank < pRank0
         % Try to increase rank up to initial rank
         [B, Q, R, E, pRank, DelVar, NewVar] = FullRank(B, A, Q, R, E, pRank,...
          pRank0, epsRank);
         if ~isempty(DelVar)
            x(DelVar)=0;
         end
         iB = find(B==0);
         iN = find(B==1);
         if isempty(iN)
            x(iB) = tomsol(6, Q, R, E, pRank, b_UEQ);
         else
            x(iB) = tomsol(6, Q, R, E, pRank, b_UEQ-A(:,iN)*x(iN));
         end
      else
         E=sparse(P(:),[1:m]',ones(m,1),m,m);
      end
   end

   if wait & PriLev > 1
      disp('Press any key to continue ...')
      pause
   end
       %if mEQ > 0
       %   x(iEQ) = x_U(iEQ);
       %end
end

%if PriLev >= 1
%   fprintf('\n--- DualSolve: TOO MANY ITERATIONS --- ITER = %7.0f \n',Iter)
%end
   fprintf('\n--- DualSolve: TOO MANY ITERATIONS --- ITER = %7.0f \n',Iter)

ResultDLP.x_k=xScale(x, Ascale, D2, n0);
ResultDLP.y_k=yScale(y, Ascale, D1, m0);
ResultDLP.f_k=f_k;
B(iEQ)=0;
ResultDLP.QP.B=B(1:n0);
ResultDLP.ExitFlag=1;
ResultDLP.ExitText=ExitText(1);

% ------------------------------------------------------------------------
function [B, Q, R, E, pRank, DelVar, NewVar] = FullRank(B, A, Q, R, E, ...
          pRank, pMax, epsRank) 
% ------------------------------------------------------------------------

nargin;

[m n]= size(A);

NewVar=[];
DelVar=[];
pRank0=pRank;

if pRank < pMax
   % Too few variables in basis or rank deficiency in linear constraints

   maxR=max(abs(diag(R)));
   [P,jE]=find(E);
   iN  = find(B == 1 | B==-1);

   if pRank==1
      iz=iN;
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
         % Save variable number for later check
         NewVar=[NewVar;k];

         if length(P) > m
            DelCol=P(m+1);
            [Q, R, P, pRank] = DeleteQR(Q, R, P, pRank, DelCol);
            k=iB(DelCol);
            % Put deleted variable on lower bound. Later check if OK.
            B(k)=1;
            % Save variable number for later check
            DelVar=[DelVar;k];

         end
      end
   end

if any(sort(P)~=(1:length(P))') 
   disp('DualSolve (FullRank): Serious P error');
   disp('Bad input data, or some logical program error.');
   size(R)
   size(E)
   length(P)
   xprinti(P,'P:');
   pause
end

   if pRank > pRank0
      % Changes has been made to the QR decomposition
      E = sparse(P,1:length(P),ones(length(P),1),length(P),length(P));
      Q = Q(1:m,1:m);
      R = R(1:m,1:pRank);
   end
end

% ------------------------------------------------------------------------
function [x] = xScale(x, Ascale, D2, n0)
% ------------------------------------------------------------------------
if Ascale
   x=x(1:n0).*D2(1:n0);
else
   x=x(1:n0);
end
% ------------------------------------------------------------------------
function [y] = yScale(y, Ascale, D1, m0)
% ------------------------------------------------------------------------
if Ascale
   y = y(1:m0)./D1;
else
   y = y(1:m0);
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
     Text = 'Infeasible dual problem';
   case 4
     Text = 'Illegal step length due to numerical difficulties';
   case 6
     Text = 'No dual feasible starting point found';
   case 7
     Text = str2mat('Illegal step length due to numerical difficulties' ...
                   ,'Should not occur!');
   case 8
     Text = 'Convergence because f_k >= QP.DualLimit';
   case 9
     Text = str2mat('x_L(i) > x_U(i) + xTol for some i' ... 
                   ,'No solution exists!');
end

% MODIFICATION LOG:
%
% 980717  mbk  'A=[A,[zeros(meq,mi),eye(mi,mi)]]' changed to
%              'A=[A,[zeros(meq,mi);eye(mi,mi)]]' on line 92.
%
% 981111  hkh  Changed Dantzig update from
%                      iN(pHat) = p; iB(qHat) = q; to
%                      iN(qHat) = p; iB(pHat) = q;
% 981111  hkh  Fix y_0 and otherwise x_0 to be used, if B not given.
% 981117  hkh  Revise for new Result/Prob input/output format.
%              Use cpTransf and QR and SVD factorization
% 981119  hkh  Change field to Prob.QP.B
% 981123  hkh  Clean up, change comments.  Set default xTol to 100*eps.
% 990810  hkh  Revised for TOMLAB v2.0. Use QR with update.
% 990811  hkh  Checking feasibility, revising print levels.
% 000916  hkh  Adding text for ExitText
% 040602  med  Help fix PriLev to PriLevOpt
% 040828  ango Check that Prob.QP.B is of correct length.
% 041222  med  Safeguard added for x_0
% 050604  hkh  Change call to lpSolve to new name lpSimplex
% 051215  med  case 5 switched to case 4 in ExitText, new text
