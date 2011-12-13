% qpSolve.m
%
% function Result = qpSolve(Prob)
%
% Active set strategi for Quadratic Programming.
%
% Minimization problem:
%
%        min   0.5 * x' * F * x + c' * x.  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%
% Equality equations: Set b_L==b_U
% Fixed    variables: Set x_L==x_U
%
% INPUT PARAMETERS
% Fields in Prob:
%   QP.F:   The matrix F in 0.5 x' F x
%   QP.c:   The vector c in c'x
%   A:      The linear constraint matrix
%   b_L:    The lower bounds for the linear constraints
%   b_U:    The upper bounds for the linear constraints
%   x_L:    Lower bounds on x
%   x_U:    Upper bounds on x
%           b_L, b_U, x_L, x_U must either be empty or of full length
%   x_0:    Starting point x
%PriLevOpt  Print level: 0 None, 1 Final result, 2 Each iteration, short
%           3 More information each iteration
%   SolverFP Name of the solver used for FP (feasible point) subproblems.
%            If empty, picked from a list, best available with a license
%
% optParam struct:
%   wait:   Pause at each iteration if wait is true ( = 1)
%   MaxIter max([500,3*n,optParam.MaxIter]);  Maximal number of iterations
%   epsRank Rank tolerance
%   xTol    Tolerance to judge that x are close
%   bTol    Linear feasibility tolerance
% -----------------------------------------------------------------------
% Hot basis (Warm start) variables:
% Prob.QP.UseHot   If > 0, Uses a warm start basis, either from file or from
%                    the structure Prob
% Prob.QP.HotFile If nonempty, a warm start basis are read from this file
%
% If Prob.QP.HotFile is empty, read warm start basis from structure Prob:
% x = Prob.QP.Hot.x; B = Prob.QP.Hot.B;
% Q = Prob.QP.Hot.Q; R = Prob.QP.Hot.R; E = Prob.QP.Hot.E;
%
% Prob.QP.HotFreq  If > 0, a warm start basis are saved, either on a file
% HotP1xxx or HotP2xxx, where HotP1 is used when solving PhaseI problems,
% and HotP2 when solving PhaseII problems. xxx is a number.
%
% Prob.QP.HotN     The number of warm start basis files, xxx=mod(Freq,HotN),
%                    where xxx is a number from 0 to HotN-1
% -----------------------------------------------------------------------
%
%
% OUTPUT PARAMETERS
% Structure Result. Fields used:
%   Iter     Number of iterations
%   ExitFlag Exit flag
%            == 0  => OK
%            == 1  => Maximal number of iterations reached. No bfs found.
%            == 2  => Unbounded feasible region.
%            == 3  => Rank problems
%            == 4  => No feasible point found with lpSimplex
%            == 10 => Errors in input parameters
%   ExitTest Text string giving ExitFlag and Inform information
%   Inform   If ExitFlag > 0, Inform=ExitFlag, otherwise Inform show type
%            of convergence:
%            0 = Unconstrained solution
%            1 = lambda >= 0.
%            2 = lambda >= 0. No 2nd order Lagrange mult. estimate available
%            3 = lambda and 2nd order Lagrange mult. positive, problem is
%                not negative definite.
%            4 = Negative definite problem. 2nd order Lagrange mult. positive,
%                but releasing variables leads to same working set.
%   x_k      Solution
%   v_k      Lagrange parameters; lower/upper bounds, then linear constraints
%   p_dx     Search steps in x
%   alphaV   Step lengths for each search step
%   f_k      Function value 0.5*x'*F*x+c'*x
%   g_k      Gradient F*x+c
%   H_k      Hessian F (constant)
%   Solver   qpSolve
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1994-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Apr 20, 1994.   Last modified Jul 17, 2009.

function ResultQP = qpSolve(Prob)

if nargin < 1
   error('qpSolve needs input structure Prob');
end

solvType=checkType('qp');
Prob=ProbCheck(Prob,'qpSolve',solvType);
Prob = iniSolveMini(Prob);

xEMPTY=isempty(Prob.x_0);
if ~xEMPTY
   xEMPTY=all(Prob.x_0==0);
end

% Define QP matrices and vectors from Prob structure

[n, x, x_km1, xEqual, x_L, x_U, Prob] = BoundInit(Prob);

xFree=isinf(x_L) & isinf(x_U);
iFree=find(xFree);

% Linear constraints

[mA, Ax, bEqual, b_L, b_U, A] = LinearConstr(Prob);

F      = Prob.QP.F;
c      = Prob.QP.c(:);
n      = max([size(F,1),length(c),size(A,2)]);
N      = mA+n;
lambda = zeros(N,1);
if isfield(Prob.QP,'B')
   B   = Prob.QP.B;
else
   B   = [];
end
if isempty(c), c=zeros(n,1); end

%m      = mA;

optParam=Prob.optParam;

if isfield(Prob.QP,'UseHot')
   if isempty(Prob.QP.UseHot),  Prob.QP.UseHot=0; end
   if isempty(Prob.QP.HotFreq), Prob.QP.HotFreq=0; end
   if isempty(Prob.QP.HotN),    Prob.QP.HotN=0; end

   UseHot  = Prob.QP.UseHot;
   HotFreq = Prob.QP.HotFreq;
   HotN    = Prob.QP.HotN;
else
   UseHot  = 0;
   HotFreq = 0;
   HotN    = 0;
end

Freq=0;

if HotFreq > 0
   HotName='HotQP';
end

PriLev    = Prob.PriLevOpt;       % Print level
IterPrint = optParam.IterPrint;   % Print short information each iteration
wait      = optParam.wait;
MaxIter   = max([500,3*n,optParam.MaxIter]);  % Maximal number of iterations

NOSTOP=1;

if isempty(Prob.SolverFP)
   if checkLicense('minos')
      SolverFP='lp-minos';
   elseif checkLicense('bqpd')
      SolverFP='bqpd';
   elseif checkLicense('xa')
      SolverFP='xa';
   else
      SolverFP='milpSolve';
   end
   % SolverFP=GetSolver('fp',Prob.LargeScale);
   Prob.SolverFP = SolverFP;
else
   SolverFP=Prob.SolverFP;
end

% Setup Structure used in QP calls
ProbQP = CreateProbQP(Prob, 'fp', max(1000,3*n), PriLev-3, optParam.wait);


%DEBUG=0;

epsRank = optParam.eps_Rank;
xTol=optParam.xTol;
bTol=optParam.bTol;


if isempty(bTol), bTol   = (1000+N+n)*eps;   end
if isempty(bTol), xTol   = 1E-12;            end

lTol=[xTol*ones(n,1);bTol*ones(mA,1)];

ResultQP                 = ResultDef(Prob);
ResultQP.Prob            = Prob;
ResultQP.Solver          = 'qpSolve';
ResultQP.SolverAlgorithm = 'Active set strategy';
ResultQP.f_0             = [];
ResultQP.H_k             = F;

if isempty(F)
   LP=1;
else
   LP=(0==sum(sum(F~=0)));
end

if LP
   F_eig=zeros(n,1);
   F_0=n;
   F_rank=0;
else
   F_rank=n;   % Rank of F. Assumed full rank
   F_0=-1;     % Number of zero eigenvalues for F. -1 == not computed
end

% Equalities in set EQ, always active. iE is index pointers to EQ
EQ=cat(1,xEqual, bEqual);

iE= find(EQ);

bl=cat(1,x_L,b_L);
bu=cat(1,x_U,b_U);


if issparse(A)
   D=[speye(n,n);A];
elseif isempty(A)
   D=speye(n,n);
else
   D=cat(1,eye(n),A);
end

rel_idx=[];
rel_B=[];
F_Q=[];F_R=[];F_E=[];F_v=[];F_eig=[];
FP=0;
Ax=zeros(0,1);
if PriLev > 2
   xprinte(c,'c:');
   format compact
   disp('F:');
   disp(F);
   rank(F)
   xprinte(eig(F),'eig:')
   if wait, pause; end
end

if UseHot
   if PriLev > 2, disp('Use Hot Basis'); end
   if ~isempty(Prob.QP.HotFile)
      load(Prob.QP.HotFile,'x','B','Q','R','E');
      if PriLev > 2
         fprintf('qpSolve: Load Hot Basis File %s',Prob.QP.HotFile)
         fprintf('\n')
      end
   else
      if PriLev > 2
         fprintf('qpSolve: Load Hot Basis From Prob structure')
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
   x = x(1:n);
   B = B(1:n);
   if mA > 0
      Ax=A*x;
   end
   g=cat(1,x,Ax);
elseif xEMPTY & ~isempty(B)
   x =zeros(n,1);
   iL=find(B(1:n)== 1);
   iU=find(B(1:n)==-1);
   iB=find(B(1:n)== 0);
   x(iL)=x_L(iL);
   x(iU)=x_U(iU);
   if mA > 0
      if length(B) < N
         B=[B(1:n);zeros(mA,1)];
         Ax=A*x;
         % Active on lower bound or in equality set EQ
         %iL= Ax-bTol <= b_L | bEqual;
         iL= Ax-b_L <= bTol*max(1,abs(Ax)) | bEqual;
         % Active on upper bound
         %iU= Ax+bTol >= b_U & ~bEqual;
         iU= b_U-Ax <= bTol*max(1,abs(Ax)) & ~bEqual;
      end
      B(n+find(iL))=1;
      B(n+find(iU))=-1;
      iAU=find(B(n+1:N)== 1);
      iAL=find(B(n+1:N)==-1);
      if isempty(iB)
         Q=[]; R=[]; E=[]; pRank=0;
      else
         % Compute the QR-decomposition
         [Q, R, E, pRank] = ComputeQR([A(iAU,iB);A(iAL,iB)], epsRank);
         if isempty(iAU)
            x(iB) = tomsol(6, Q, R, E, pRank, b_L(iAL)-A(iAL,:)*x);
         elseif isempty(iAL)
            x(iB) = tomsol(6, Q, R, E, pRank, b_U(iAU)-A(iAU,:)*x);
         else
            x(iB) = tomsol(6, Q, R, E, pRank, ...
                        [b_U(iAU)-A(iAU,:)*x; b_L(iAL)-A(iAL,:)*x]);
         end
      end
   end
   g=cat(1,x,Ax);

   if PriLev > 1
      fprintf('Computed starting value x as\n');
      xprint(x,'x:')
   end
else
   iEq=find(bEqual);
   if ~isempty(iFree) & ~isempty(iEq) & xEMPTY
      % Try to estimate a better x solution
      x(iFree)=0;
      x(iFree)=A(iEq,iFree) \ (b_U(iEq) - A(iEq,:)*x);
      if PriLev > 1
         disp('Free var solution')
         xprinte(x(iFree),'x:');
      end
   end
   if mA > 0
      Ax=A*x;
   end
   g=cat(1,x,Ax);
   % Active on lower bound or in equality set EQ
   %iL= g-lTol <= bl | EQ;
   iL = g  - bl <= lTol .* max(1,abs(g)) | EQ;
   % Active on upper bound
   %iU= g+lTol >= bu & ~EQ;
   iU = bu - g  <= lTol .* max(1,abs(g)) & ~EQ;

   B=zeros(N,1);
   B(find(iL))=1;
   B(find(iU))=-1;
   if (~isempty(iFree) & length(iFree) < n) & isempty(iEq)
      FP=1;
   end
end

p_V = [];
alfa_V = [];

% Check if x is feasible

%if DEBUG
%   %disp('bl g bu')
%   xprint(x,'New x: ');
%   PrintMatrix([bl g bu],'bl x;Ax bu',120)
%   if wait, pause; end
%end

NonFeas = any( Ax-b_L < -bTol*max(1,abs(Ax)) | ...
               b_U-Ax < -bTol*max(1,abs(Ax)));

if NonFeas | FP
   if NonFeas & PriLev > 1
      disp('qpSolve: initial x not feasible')
   end
   if NonFeas & PriLev > 2
      for i=1:mA
          fprintf('%d: b_L%20.14f Ax%20.14f b_U%20.14f\n',...
                   i,b_L(i),Ax(i),b_U(i));
      end
   end

   x_0      = x;

   ResultLP = DoFP(SolverFP,ProbQP,Prob);

   ExitFlag = ResultLP.ExitFlag;
   x        = ResultLP.x_k;

   if ExitFlag > 0  %  Phase 1 failed
      ResultQP.Iter     = 0;
      ResultQP.ExitFlag = 4;
      ResultQP.Inform   = ResultLP.ExitFlag;
      ResultQP.ExitText = ExitText(4,ResultQP.Inform);
      ResultQP.f_k      = NaN;
      ResultQP.g_k      = [];
      ResultQP.x_k      = x;
      ResultQP.QP.B     = ResultLP.QP.B;
      ResultQP.r_k      = ResultLP.r_k;
      ResultQP.y_k      = ResultLP.y_k;
      ResultQP.v_k      = lambda;
      % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
      if ~isempty(x)
         % ResultQP.xState=double( (x==x_L)+2*(x==x_U) );
         ResultQP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );
      end
      if isfield(ResultLP.QP,'Hot');
      end
      ResultQP = endSolveMini(Prob,ResultQP);
      return
   end

   v_k      = ResultLP.v_k;
   B        = ResultLP.QP.B(1:n);
   f_k      = ResultLP.f_k;
%f_k
   if mA > 0
      Ax=A*x;
      if length(B) < N
         B=[B(1:n);zeros(mA,1)];
         % Active on lower bound or in equality set EQ
         %iL= Ax-bTol <= b_L | bEqual;
         iL= Ax-b_L <= bTol*max(1,abs(Ax)) | bEqual;
         % Active on upper bound
         %iU= Ax+bTol >= b_U & ~bEqual;
         iU= b_U-Ax <= bTol*max(1,abs(Ax)) & ~bEqual;
      end
%Ax
%[b_L Ax b_U]
%pause
      B(find(iL))=1;
      B(find(iU))=-1;
%B
   end
   g=cat(1,x,Ax);
%g
   % Active on lower bound or in equality set EQ
   %iL= g-lTol <= bl | EQ;
   iL = g  - bl <= lTol .* max(1,abs(g)) | EQ;
   % Active on upper bound
   %iU= g+lTol >= bu & ~EQ;
   iU = bu - g  <= lTol .* max(1,abs(g)) & ~EQ;

   B=zeros(N,1);

   B(find(iL))=1;

   B(find(iU))=-1;

   p_V = x-x_0; % Reset p_V and alfa_V
   alfa_V = 1;
end

if any( Ax - b_L < -bTol*max(1,abs(Ax)))| ...
   any( Ax - b_U >  bTol*max(1,abs(Ax)))| ...
   any(  x - x_L < -xTol*max(1,abs(x))) | ...
   any(  x - x_U >  xTol*max(1,abs(x)))


%if DEBUG
%   xprinte(b_L,'b_L')
%   xprinte(b_U,'b_U')
%   xprinte(Ax,'Ax')
%   pause
%z1= (Ax - b_L) ./ max(1,abs(Ax));
%z2= (Ax - b_U) ./ max(1,abs(Ax));
%z3= (x - x_L)  ./ max(1,abs(x));
%z4= (x - x_U)  ./ max(1,abs(x));
%format long
%[mz1, iz1]=min(z1);
%[mz2, iz2]=max(z2);
%mz1
%mz2
%fprintf('%d: b_L%20.14f Ax%20.14f b_U%20.14f\n',...
%         iz1,b_L(iz1),Ax(iz1),b_U(iz1));
%fprintf('%d: b_L%20.14f Ax%20.14f b_U%20.14f\n',...
%        iz2,b_L(iz2),Ax(iz2),b_U(iz2));
%min(z3)
%max(z4)
%plot(z1)
%pause
%plot(z2)
%pause
%plot(z3)
%pause
%plot(z4)
%pause
%end

if NOSTOP
   if PriLev > 1
      disp('ACCURACY PROBLEM in PHASE1, but continue');
   end
else
   if PriLev > 1
      disp('Initial x not feasible w.r.t. bounds. No solution exists ')
   end
   if PriLev > 1
      for i=1:mA
          fprintf('%d: b_L%20.14f Ax%20.14f b_U%20.14f\n',...
                   i,b_L(i),Ax(i),b_U(i));
      end
      if wait, pause; end
      for i=1:n
          fprintf('%d: x_L%20.14f Ax%20.14f x_U%20.14f\n',...
                   i,x_L(i),x(i),x_U(i));
      end
      if wait, pause; end
   end
   ResultQP.Iter=0;
   ResultQP.ExitFlag=4;
   ResultQP.Inform=5;
   ResultQP.ExitText = ExitText(4,5);
   ResultQP.f_k=NaN;
   ResultQP.g_k=[];
   ResultQP.x_k=x;
   ResultQP.QP.B=[];
   ResultQP.v_k = lambda;
   % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
   %ResultQP.xState=double( (x==x_L)+2*(x==x_U) );
   ResultQP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );
   ResultQP = endSolveMini(Prob,ResultQP);
   return
end
end

if isempty(F)
   f_k=c'*x;
   g_k=c;
else
   f_k=0.5*(x'*F*x)+c'*x;
   g_k=F*x+c;
end

ResultQP.x_k=x;
ResultQP.f_k=f_k;
ResultQP.g_k=g_k;

ResultQP.x_0=x;
ResultQP.f_0=f_k;

% Free==0; On lower == 1; On upper == 2; Fixed == 3;
%ResultQP.xState=double( (x==x_L)+2*(x==x_U) );

ResultQP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );

k = 1;
jAdd=[];
P=[];

while 1
   % Determine working set of active constraints

   W=abs(B) | EQ;

   nAct=sum(W);

%if DEBUG
%   z=min(g-bl,bu-g);
%   xprinte(z(find(W)),'W: ');
%   xprinte(z(find(~W)),'~W:');
%end

   if PriLev > 2, fprintf('No of active constraints %d\n',nAct); end

   % Flag if eigenvalue direction used
   eig_p=0;

   % Flag if gradient step
   gradStep=0;

   % Flag if negative definite problem
   NEGDEF=0;

   % Gradient in current point x
   if isempty(F)
      g0 = c;
   else
      g0 = F*x+c;
   end
   if IterPrint
      %f_k=0.5*x'*F*x+c'*x;
      fprintf('Iter %d. f(x) %20.15f\n',k,f_k);
   end

   %alphaCon=1;

   if nAct==0                         % Unconstrained solution
      [F_R F_p]=chol(F);

      if F_p > 0 % Negative definite problem
         if PriLev > 2, disp('F Negative Definite.');end

         % Check if eigenvalus to F already computed
         if size(F_eig,1)==0
            [F_v,F_eig]=eig(full(F)); % Compute eigenvalues(eigenvectors) of F
            F_eig = real(F_eig);
            F_v   = real(F_v);
            F_0   = sum(abs(diag(F_eig)) <= epsRank);
            F_rank=n-F_0;
            if PriLev > 3
               fprintf('F rank: %4.0f. (Empty Working Set)\n',F_rank);
            end
         end
         [F_min_eig jE]=sort(diag(F_eig));
         if F_min_eig(1) < -bTol, NEGDEF=1;end
         %p=F_v(:,jE(1)); % Find direction of highest negative curvature

         % Try up to five computed directions!
         p=F_v(:,jE(min(5,length(jE))));

         p_g0=p'*g0; % Directed derivative in eigenvalue directions

         ip = find(p_g0 > 0);

         if ~isempty(ip)
            p(:,ip)=-p(:,ip);    % Change sign to get descent
            p_g0(:,ip)=-p_g0(:,ip);
         end
         if norm(g0) > xTol
            p=[p, -g0];
         end

         eig_p=1;
      else
         % Use safe projection solution in case of rank defiency
         if isempty(F)
            p = -g0;
            gradStep=1;
         else
            zz = diag(F_R);
            condF=min(zz)/max(zz);
            if condF < 1E-15
%condF
               % F is rank deficient
               [Q1, R1, E1, pRank] = ComputeQR(F,epsRank);
               p = tomsol(6, Q1, R1, E1, pRank, -g0);

            else
               p = F_R \(F_R'\(-g0));
            end
         end
      end
      p_sq=max(1E-12,p'*p); % Convergence now. p_sq >= 1E-12 ==> Conv.test

   else  % There are constraints present
      %
      % Solve QPE:     min  0.5*p'*F*p+c'*p
      %                s/t  D*p  = 0
      %

      % Build constraint matrix and find rank pRank

      % H=diag(B(find(B~=0)))*D(find(W),:);
      [Q, R, E, pRank] = ComputeQR([diag(B(find(B~=0)))*D(find(W),:)]',epsRank);

      if pRank < n             % Not full rank
         Z = Q(:,pRank+1:n);   % An orthogonal base for the null space

         Zg = Z'* g0;          % Projected gradient

         if ~isempty(F)
            ZFZ = Z'*F*Z;      % Projected Hessian
            [QZ, RZ, EZ, pRankZ] = ComputeQR(full(ZFZ),epsRank);

            [Z_R Z_p]=chol(ZFZ);

            %if size(ZFZ,1)==1 &  abs(ZFZ(1,1)) < 1E-12
            %   % Too small value to be considered full rank %   Z_p=NaN;
            %end

         else
            Z_p=0;
         end

         if isnan(Z_p)         % Too small value
            gradStep=1;
            p = -Z*Zg;
            %p = zeros(n,1);
         elseif Z_p > 0            % Negative definite
            [V,L]=eig(full(ZFZ));
            L = real(L);
            V = real(V);
            [min_eig jE]=sort(diag(L));
            if min_eig(1) < -bTol, NEGDEF=1;end
            if PriLev > 2, disp('Negative curvature'); end

            xsi=V(:,jE(1));
            %p=Z*xsi; % Find direction of negative curvature

            % Try up to five computed directions!
            p=Z*V(:,jE(min(5,length(jE))));

            p_g0=p'*g0; % Directed derivative in eigenvalue direction

            ip = find(p_g0 > 0);

            if ~isempty(ip)
               p(:,ip)=-p(:,ip);    % Change sign to get descent
               xsi=-xsi;
               p_g0(:,ip)=-p_g0(:,ip);
            end

            if norm(Zg) > xTol
               p=[p, -Z*Zg];
            end

            %p=[p, -Z*Zg, -Z*Z'*eye(n), eye(n)];

            % Scale the step

            XSIZFZ=xsi'*ZFZ*xsi;
            if XSIZFZ > xTol
               tau=-xsi'*Zg/XSIZFZ;
               p=tau*p;
            end

            eig_p=1;

         else
            % Solve reduced linear system
            % Null space part is solution x
            if isempty(F)
               gradStep=1;
               p = -Z*Zg;

%if DEBUG
%   disp('null space sol')
%   p
%end

            else
               zz = diag(Z_R);
               condZFZ=min(zz)/max(zz);
               if condZFZ < 1E-15
%condZFZ
                  % ZFZ is rank deficient
                  [Q1, R1, E1, pRank] = ComputeQR(ZFZ,epsRank);
                  p = Z * tomsol(6, Q1, R1, E1, pRank, -Zg);

               else
                  p = Z * (Z_R \(Z_R'\(-Zg)));
               end

%if DEBUG
%   disp('Z_R solution')
%   p
%   norm(p)
%end



               if norm(p) > 1E5*max(1,norm(x))
                  gradStep=1;
                  p = -Z*Zg;
               end

            end
         end
         if size(p,2) > 1
            p_sq=1;
         else
            p_sq=max(p'*p);
         end
      elseif LP
         p=[];
         p_sq=0;
      else % Full rank ==> p==0
         if F_0 < 0
            [F_v F_eig]=eig(full(F));
            F_eig = real(F_eig);
            F_v   = real(F_v);
            F_0   = sum(abs(diag(F_eig)) <= epsRank);
            F_rank= n-F_0;
            if PriLev > 3
               fprintf('F rank: %4.0f. (Full rank in A) \n',F_rank);
            end
         end
         [F_min_eig j]=min(diag(F_eig));
         if F_min_eig < -bTol, NEGDEF=1;end
         p=[];
         p_sq=0;
      end
   end

   if PriLev > 1 & length(p) > 0
      for ij=1:size(p,2)
          fprintf('Search direction #%d p:\n',ij);
          xprint(p(:,ij),'p:  ');
      end
   end

   if p_sq < 1E-12
      % Estimate Lagrange multipliers

      lambda=zeros(N,1);
      if pRank > 0
         l = tomsol(6, Q, R, E, pRank, g0);
         lambda(find(W))=l;
         % HKH FIX
         if all(lambda >= 0)
            ix = find(g(n+1:end) > bu(n+1:end) & lambda(n+1:end) >= 0);
            if ~isempty(ix)
               %g(n+ix)
               %bu(n+ix)
               %lambda(n+ix)
               %ix
               lambda(n+ix) = -lambda(n+ix);
            end
            ix = find(g(n+1:end) < bl(n+1:end) & lambda(n+1:end) >= 0);
            if ~isempty(ix)
               %g(n+ix)
               %bl(n+ix)
               %lambda(n+ix)
               %ix
               lambda(n+ix) = -lambda(n+ix);
            end
         end
      end

      if PriLev > 1
         fprintf('Lagrange multipliers lambda for equalities:\n');
         xprint(lambda(iE),'l_E:  ');
         fprintf('Lagrange multipliers for inequalities:\n');
         xprinte(lambda(~EQ),'l_I:');
      end

      if PriLev > 2, disp('Check for Negative Lagrange multipliers'); end

      q=find((lambda < -lTol.*max(1,abs(g))) & (W & ~EQ));
      if PriLev > 2
         disp('Check for Negative Lagrange multipliers on x_L and x_U');
      end
      if isempty(q)
         if PriLev > 2, disp('Check for zero multipliers'); end

         % No negative multipliers, but maybe zero multipliers
         q0=find(abs(lambda) < lTol .* max(1,abs(g)) & (W & ~EQ));
         if ~isempty(q0) & pRank > 0
            if PriLev > 2
               disp('Check zero multipliers using 2nd order estimate');
            end
            % Compute 2nd order multiplier estimates, lambda2
            q2=[];
            if size(p_V,2) > 0 & size(F,1)==n
               % Use the previous search direction, last column in p_V

               eta = tomsol(6, Q, R, E, pRank, g0+F*p_V(:,size(p_V,2)));

               %if issparse(H)
               %   eta=(R(1:pRank,1:pRank)\...
               %      (Q(:,1:pRank)'*(g0+F*p_V(:,size(p_V,2)))));
               %else
               %   eta=P(:,1:pRank)*(R(1:pRank,1:pRank)\...
               %      (Q(:,1:pRank)'*(g0+F*p_V(:,size(p_V,2)))));
               %end

               lambda2=zeros(N,1);
               lambda2(find(W))=eta;
               %q2=find(lambda2 < -lTol & (W & ~EQ));
               q2=find(lambda2 < -lTol.*max(1,abs(g)) & (W & ~EQ));
               if isempty(q2) & ~NEGDEF
                  % Convergence, both 1st and 2nd order Lagrange multiplier
                  % estimates are non negative
                  if PriLev >= 1,disp('lambda and lambda2 >=0'),end
                  ResultQP.Iter=k;
                  ResultQP.ExitFlag=0;
                  ResultQP.Inform=1;
                  ResultQP.ExitText = ExitText(0,1);
                  ResultQP.QP.B=B(1:n);
                  ResultQP.v_k = lambda;
                  if mod(k,HotFreq)~=0 & HotFreq > 0
                     Freq=Freq+1;
                     save([HotName num2str(mod(Freq,HotN))],...
                         'x','B','Q','R','P');
                  end
                  ResultQP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );
                  ResultQP = endSolveMini(Prob,ResultQP);
                  return;
               elseif ~isempty(q2)
                  if PriLev >= 2,disp('Use lambda2 to delete'),end
                  % Delete from working set
                  %if 1
                     B(q2)=0;
                  %else
                  %   % Cautious strategy, only delete one constraint
                  %   [q_min, qix]=min(lambda2(q2))
                  %   B(q2(qix))=0;
                  %end % if 0/1
               end
            end
            if NEGDEF & isempty(q2)
               % Try releasing all zero multiplier variables
               if length(q0)==length(rel_idx)
                  if all(q0==rel_idx)
                     if PriLev >= 1,disp('lambda >=0. Same working set'),end
                     ResultQP.Iter=k;
                     ResultQP.ExitFlag=0;
                     ResultQP.Inform=2;
                     ResultQP.ExitText = ExitText(0,2);
                     ResultQP.QP.B=B(1:n);
                     ResultQP.v_k = lambda;
                     if mod(k,HotFreq)~=0 & HotFreq > 0
                        Freq=Freq+1;
                        save([HotName num2str(mod(Freq,HotN))],...
                            'x','B','Q','R','P');
                     end
                     ResultQP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );
                     ResultQP = endSolveMini(Prob,ResultQP);
                     return;
                  end
               end
               if PriLev > 2,disp('Releasing all zero Lagrange multipliers');end

               if ~isempty(q0)
                  B(q0)=0;
                  rel_idx=q0;   % Save variable index for later test
               end
            else
               if PriLev >= 1,disp('lambda >=0. No estimate on lambda2'),end
               ResultQP.Iter=k;
               ResultQP.ExitFlag=0;
               ResultQP.Inform=3;
               ResultQP.ExitText = ExitText(0,3);
               ResultQP.QP.B=B(1:n);
               ResultQP.v_k = lambda;
               if mod(k,HotFreq)~=0 & HotFreq > 0
                  Freq=Freq+1;
                  save([HotName num2str(mod(Freq,HotN))],...
                      'x','B','Q','R','P');
               end
               ResultQP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );
               ResultQP = endSolveMini(Prob,ResultQP);
               return;
            end
         else
            if PriLev >= 1,disp('lambda >=0.'),end
            ResultQP.Iter=k;
            ResultQP.ExitFlag=0;
            ResultQP.Inform=1;
            ResultQP.ExitText = ExitText(0,1);
            ResultQP.QP.B=B(1:n);
            ResultQP.v_k = lambda;
            if mod(k,HotFreq)~=0 & HotFreq > 0
               Freq=Freq+1;
               save([HotName num2str(mod(Freq,HotN))],...
                   'x','B','Q','R','P');
            end
            ResultQP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );
            ResultQP = endSolveMini(Prob,ResultQP);
            return
         end
      elseif LP
         % LP search direction
         j1=sum(W(1:q(1))~=0);
         % LP search direction
         p=(R(1:pRank,1:pRank)\Q(:,1:pRank)')'*E(j1,1:pRank)';

         %p2 = tomsol(6, Q, R, E, pRank, );
         % H=diag(B(find(B~=0)))*D(find(W),:);

         if PriLev > 1
            fprintf('Search direction p:\n');
            xprint(p,'p:  ');
         end
         % Slower to compute ... p=A_inv(:,me+q(1));
         B(q(1))=0;  % Constraint q(1) deleted from working set
      elseif sum(W)+F_rank < n+1
         % Rank problem. No deletion is possible
         %c1=g-100*lTol <= bl & ~W;
         %c2=g+100*lTol >= bu & ~W;
         c1=g-bl <= 100*lTol.*max(1,abs(g))  & ~W;
         c2=bu-g <= 100*lTol.*max(1,abs(g))  & ~W;
         if PriLev > 2
            disp('RANK PROBLEM! Check for active constraints in inactive set');
         end
         if any(c1 | c2)
            g1=Inf*ones(N,1);
            g2=g1;
            g1(find(c1))=g(find(c1))-bl(find(c1));
            [j1min j1]=min(g1);
            g2(find(c2))=bu(find(c2))-g(find(c2));
            [j2min j2]=min(g2);
            if j1min <= j2min
               B(j1)=1;   % Add constraint j1 to active set
            else
               B(j2)=-1;   % Add constraint j1 to active set
            end
            B(q(1))=0;  % Constraint q(1) deleted from working set
         else
            if PriLev >= 0
               disp('Rank problem. CAN NOT FIND SOLUTION TO QP!!!');
            end
            ResultQP.Iter=k;
            ResultQP.ExitFlag=3;
            ResultQP.Inform=3;
            ResultQP.ExitText = ExitText(3,3);
            ResultQP.QP.B=B(1:n);
            ResultQP.v_k = lambda;
            if mod(k,HotFreq)~=0 & HotFreq > 0
               Freq=Freq+1;
               save([HotName num2str(mod(Freq,HotN))],...
                   'x','B','Q','R','P');
            end
            ResultQP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );
            ResultQP = endSolveMini(Prob,ResultQP);
            return
         end
      else
         % Take 1st index with negative Lagrange multiplier
         % BUT NOT THE LAST INCLUDED, if only one was added
         if length(q)==1 | length(jAdd) ~= 1
            B(q(1))=0;  % Constraint q(1) deleted from working set
               if length(q) > 1
                  %ixx=find(lambda(q) < -1);
                  %if ~isempty(ixx)
                  %   B(q(ixx))=0;
                  %   disp('FIX Release all')
                  %end
                  % If the 2nd largest lambda is big, release also that one
                  if lambda(q(2)) < -1
                     B(q(2))=0;
                  end
               end
         else
            if q(1) == abs(jAdd)
               B(q(2))=0;  % Constraint q(2) deleted from working set
            else
               B(q(1))=0;  % Constraint q(1) deleted from working set
            end
         end
      end
   end

% 070726 frhe Modified expression to behave like before the | to |-change
   if all(p_sq(:) >= 1E-12) | LP % STEP in subspace of Working Set.

      % Previous point feasible or not
      NonFeas = any( Ax-b_L < -bTol*max(1,abs(Ax)) | ...
                     b_U-Ax < -bTol*max(1,abs(Ax)));

      [alpha, p, j] = DetermineAlpha(p, W, g0, D, lTol, F, c, g,...
                                     bl, bu, x, eig_p, LP, gradStep);
      x = x + alpha*p;

      gOld = g;

      if mA > 0, Ax=A*x; end
      g=[x;Ax];

      f_km1=f_k;
      if isempty(F)
         f_k=c'*x;
         g_k=c;
      else
         f_k=0.5*(x'*F*x)+c'*x;
         g_k=F*x+c;
      end

      ResultQP.x_k=x;
      ResultQP.f_k=f_k;
      ResultQP.g_k=g_k;
      % Free==0; On lower == 1; On upper == 2; Fixed == 3;
      %ResultQP.xState=double((x==x_L)+2*(x==x_U));
      ResultQP.xState=double((B==1) + 2*(B==-1) + 3*(B==2));

      %if f_k-f_km1 > 1E-12
      if f_k-f_km1 > 1E-12 & ~NonFeas
      if PriLev >= 1
            fprintf('WARNING: No descent in qpSolve: f_k-f_km1 =%22.15e\n',...
	             f_k-f_km1);
	     end
         x = x - alpha*p;
         g = gOld;
         alpha = 0;

         if 1
         % HKH FIX. Try quiting here. Take previous point
         ResultQP.x_k=x-alpha*p;
         ResultQP.f_k=f_km1;
         ResultQP.Iter=k;
         ResultQP.ExitFlag=0;
         ResultQP.Inform=5;
         ResultQP.ExitText = ExitText(0,5);
         ResultQP.QP.B=B(1:n);
         ResultQP.v_k = lambda;
         if mod(k,HotFreq)~=0 & HotFreq > 0
            Freq=Freq+1;
            save([HotName num2str(mod(Freq,HotN))],...
                'x','B','Q','R','P');
         end
         ResultQP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );
         ResultQP = endSolveMini(Prob,ResultQP);
         return
         else
            ResultQP.x_k=x-alpha*p;
         end
      end
      if PriLev >= 2
         fprintf('New solution x; alpha = %15.5e ',alpha);
         fprintf('Func.value = %25.16f:\n',f_k);
         xprint(x,'x:  ');
      end
      if isinf(alpha)
         if PriLev >= 1, disp('Unbounded solution to QP-problem!!!'); end
         ResultQP.Iter=k;
         ResultQP.ExitFlag=2;
         ResultQP.Inform=2;
         ResultQP.ExitText = ExitText(2,2);
         ResultQP.QP.B=B(1:n);
         ResultQP.v_k = lambda;
         if mod(k,HotFreq)~=0 & HotFreq > 0
            Freq=Freq+1;
            save([HotName num2str(mod(Freq,HotN))],...
                'x','B','Q','R','P');
         end
         ResultQP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );
         ResultQP = endSolveMini(Prob,ResultQP);
         return
      end
      p_V = [p_V,p];
      alfa_V = [alfa_V;alpha];
      jAdd=[];
      if ((alpha < 1) | eig_p | LP | gradStep) & j ~= 0
         B(abs(abs(j)))=sign(j);
         jAdd=j;
      elseif sum(W)==0
         if PriLev >= 1, disp('Unconstrained QP-optimum!!!'); end
         ResultQP.Iter=k;
         ResultQP.ExitFlag=0;
         ResultQP.Inform=0;
         ResultQP.ExitText = ExitText(0,0);
         ResultQP.QP.B=B(1:n);
         ResultQP.v_k = lambda;
         if mod(k,HotFreq)~=0 & HotFreq > 0
            Freq=Freq+1;
            save([HotName num2str(mod(Freq,HotN))],...
                'x','B','Q','R','P');
         end
         ResultQP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );
         ResultQP = endSolveMini(Prob,ResultQP);
         return
      end
      %if 1
      % Check if any new constraints have reached the bounds, and are
      % possible to add without rank problems

      W=abs(B) | EQ;

      if sum(W)+1 < n
         %HKH disp('--------------HKH-------------------')
         % It is possible to add more constraints
         c1 = g  - bl <= lTol .* max(1,abs(g)) & ~W;
         c2 = bu - g  <= lTol .* max(1,abs(g)) & ~W;
         %c1=g-lTol <= bl & ~W;
         %c2=g+lTol >= bu & ~W;
         if any(c1 | c2)
            if PriLev > 2
               disp('Adding extra constraints!!!');
            end
           %if 0
           % g1=Inf*ones(N,1);
           % g2=g1;
           % g1(find(c1))=g(find(c1))-bl(find(c1));
           % [j1min j1]=min(g1);
           % g2(find(c2))=bu(find(c2))-g(find(c2));
           % [j2min j2]=min(g2);
           % if j1min <= j2min
           %    B(j1)=1;   % Add constraint j1 to active set
           %    if PriLev > 1
           %       fprintf('Add LOW %d %d\n',j1, B(j1));
           %    end
           % else
           %    B(j2)=-1;  % Add constraint j1 to active set
           %    if PriLev > 1
           %       fprintf('Add UPP %d %d\n',j2, B(j2));
           %    end
           % end
           %else   % ADD ALL
              kk=sum(W)+1;
              i=1;
              ic1=find(c1);
              while kk<n & i < length(ic1)
                  B(ic1(i))=1;
                  if PriLev > 1
                     fprintf('Add LOW %d %d\n',ic1(i), B(ic1(i)));
                  end
                  kk=kk+1;
                  i=i+1;
                  jAdd=[jAdd;ic1(i)];
              end
              i=1;
              ic2=find(c2);
              while kk<n & i < length(ic2)
                  B(ic2(i))=1;
                  if PriLev > 1
                     fprintf('Add UPP %d %d\n',ic2(i), B(ic2(i)));
                  end
                  kk=kk+1;
                  i=i+1;
                  jAdd=[jAdd;-ic2(i)];
              end
           %end % if 0/1
         end
      end
      %end % if 0/1
   end
   k = k+1;
   if mod(k,HotFreq)==0
      Freq=Freq+1;
      save([HotName num2str(mod(Freq,HotN))],'x','B','Q','R','P');
   end
   if k > MaxIter
      ResultQP.Iter=k;
      ResultQP.ExitFlag=1;
      ResultQP.Inform=1;
      ResultQP.ExitText = ExitText(1,1);
      ResultQP.QP.B=B(1:n);
      ResultQP.v_k = lambda;
      if mod(k,HotFreq)~=0 & HotFreq > 0
         Freq=Freq+1;
         save([HotName num2str(mod(Freq,HotN))],'x','B','Q','R','P');
      end
      ResultQP.xState=double( (B==1) + 2*(B==-1) + 3*(B==2) );
      ResultQP = endSolveMini(Prob,ResultQP);
      return
   end
end

%------------------------------------------------------------------------
function ResultLP=DoFP(SolverFP,ProbQP,Prob)
%------------------------------------------------------------------------

%DEBUG=0;

%ProbQP.optParam.eps_f  = 1E-8;   % Convergence tolerance for LP

ProbQP.A    = Prob.A;
ProbQP.b_L  = Prob.b_L;
ProbQP.b_U  = Prob.b_U;
ProbQP.x_L  = Prob.x_L;
ProbQP.x_U  = Prob.x_U;
ProbQP.mLin = size(ProbQP.A,1);
ProbQP.QP.c = zeros(Prob.N,1); % For nonlinear SolverFP's

ResultLP=tomRunMini(SolverFP,ProbQP);

%if DEBUG
%   CheckFP(ResultLP,ProbQP);
%end

%------------------------------------------------------------------------
function CheckFP(ResultLP,ProbQP)
%------------------------------------------------------------------------

ExitFlag = ResultLP.ExitFlag;

x1        = ResultLP.x_k;

% 1 = QPOPT, 2 = MINOS, 3 = lpSimplex

n= length(x1);

if 1
% Run QPOPT to test results against lpSimplex on FP
Solver2 = 'qpopt';
% Run MINOS to test results against lpSimplex on FP
% Solver2 = 'lp-minos';

ResultLP2=tomRunMini(Solver2,ProbQP);

ExitFlag2= ResultLP2.ExitFlag;

x2       = ResultLP2.x_k;

if any(abs(x1-x2) > 100*eps) & (ExitFlag == 0 | ExitFlag2==0)
   fprintf('FP: diff between %s and lpSimplex\n',Solver2);
   xprinte(x1,'x1: ');
   xprinte(x2,'x2: ');
   fprintf('f1 %f f2 %f\n',ResultLP.f_k,ResultLP2.f_k);
   xprinte(ProbQP.x_L,'x_L:');
   xprinte(ProbQP.x_U,'x_U:');
   xprinte(x1-x2,'Diff:');
   Ax1=ProbQP.A*x1;
   Ax2=ProbQP.A*x2;
   e1=1E-10;

   r1L=Ax1-ProbQP.b_L;
   iz1=r1L >= -e1*max(1,ProbQP.b_L);
   r1U=Ax1-ProbQP.b_U;
   iz1=iz1 & r1U <= e1*max(1,ProbQP.b_U);

   r2L=Ax2-ProbQP.b_L;
   iz2=r2L >= -e1*max(1,ProbQP.b_L);
   r2U=Ax2-ProbQP.b_U;
   iz2=iz2 & r2U <= e1*max(1,ProbQP.b_U);

   Feas1=all(iz1);
   Feas2=all(iz2);
   Inform2=ResultLP2.Inform;
   fprintf('Exitflags %d - %d ',ExitFlag, ExitFlag2);
   fprintf('Inform2 %d  Feas1 %d Feas2 %d\n',...
           Inform2, Feas1,Feas2);
   if ~Feas1 | ~Feas2
      disp('******************* NOT FEASIBLE ********************')
      [ProbQP.b_L Ax1 Ax2 ProbQP.b_U]
      ix1=find(~iz1);
      ix2=find(~iz2);
      if ~isempty(ix1)
         for i=1:length(ix1)
             fprintf('%20.15e ',[ProbQP.b_L(ix1(i)) Ax1(ix1(i))...
                                 Ax2(ix1(i)) ProbQP.b_U(ix1(i))]);
             fprintf('\n');
         end
      end
      if ~isempty(ix2)
         for i=1:length(ix2)
             fprintf('%20.15e ',[ProbQP.b_L(ix2(i)) Ax1(ix2(i))...
                                 Ax2(ix2(i)) ProbQP.b_U(ix2(i))]);
             fprintf('\n');
         end
      end
      pause
   end
end

end

%------------------------------------------------------------------------
function [alpha, p, j] = DetermineAlpha(p, W, g0, D, lTol, F, c, g,...
                                        bl, bu, x, eig_p, LP, gradStep)
%------------------------------------------------------------------------
pp=p;
%I = find(~W|W);
I=1:length(W);

for i=1:size(pp,2)
    if size(pp,2) > 1
       p=pp(:,i);
       if i==1
          alphaTry=[];
          jSave=[];
          fTry=[];
       end
       if g0'*p > 0
          p=-p;
          pp(:,i)=-pp(:,i);
       end
    end
    if eig_p | LP | gradStep
       alpha = Inf;
    else
       alpha = 1;
       %alpha = alphaCon;
    end

    if isempty(D)
       j = 0;
       s = [];
    else
       s = D(I,:)*p;
       Step = Inf*ones(length(I),1);
       L=find(~isinf(bl(I)) & s < -lTol(I) .* max(1,abs(g(I))));
       U=find(~isinf(bu(I)) & s >  lTol(I) .* max(1,abs(g(I))));
       if ~isempty(L)
          Step(L)=(bl(I(L))-g(I(L)))./s(L);
       end
       if ~isempty(U)
          Step(U)=(bu(I(U))-g(I(U)))./s(U);
       end
       [alpha0 j]=min(Step);
       if ~isempty(U)
          if any(j==U), j=-j; end
       end
       if alpha <= alpha0
          j=0;
       else
          alpha=alpha0;
       end
    end
    if size(pp,2) > 1
       alphaTry=[alphaTry;alpha];
       xTry=x+alpha*p;
       jSave=[jSave;j];
       if isempty(F)
          f_k=c'*xTry;
       else
          f_k=0.5*(xTry'*F*xTry)+c'*xTry;
       end
       fTry=[fTry;f_k];
    end
end

if size(pp,2) > 1

%fTry
%alphaTry
%jSave

   [minF j1]=min(fTry);
   %ix=find(minF==fTry);
   %if length(ix) > 1
   %   [mina ixa]=sort(alphaTry(ix));
   %   j1=ix(ixa(length(ix)));
   %end
   p=pp(:,j1);
   alpha=alphaTry(j1);
   j=jSave(j1);
end
% if I not all constraints, then j=I(j), but with correct sign.
%j1
%j
%pause

%   ExitFlag Exit flag
%            == 0  => OK
%            == 1  => Maximal number of iterations reached. No bfs found.
%   Inform   If ExitFlag > 0, Inform=ExitFlag, otherwise Inform show type
%            of convergence:
%            0 = Unconstrained solution
%            1 = lambda >= 0.
%            2 = lambda >= 0. No 2nd order Lagrange mult. estimate available
%            3 = lambda and 2nd order Lagrange mult. positive, problem is
%                not negative definite.
%            4 = Negative definite problem. 2nd order Lagrange mult. positive,
%                but releasing variables leads to same working set.
% ------------------------------
function Text = ExitText(ExitFlag,Inform)
% ------------------------------

switch  ExitFlag
   case 0
     switch  Inform
        case 0
          Text = 'Unconstrained solution';
        case 1
          Text = 'First order multipliers >= 0';
        case 2
          Text = 'First order multipliers >= 0. No 2nd order estimate used';
        case 3
          Text = 'First and 2nd order multipliers positive';
        case 4
          Text = 'Negative definite problem. Cannot proceed more';
        case 5
          Text = 'No descent. Cannot proceed more';
     end
     Text = str2mat('Optimal point found',Text);
   case 1
     Text = 'Maximal number of iterations reached';
   case 2
     Text = 'Unbounded feasible region';
   case 3
     Text = 'Rank problems';
   case 4
     Text = 'No feasible point found with lpSimplex';
   case 10
     Text = 'Errors in input parameters';
end


% MODIFICATION LOG:
%
% 980825  hkh  NOT OK: tau=-xsi'*Z'*Zg/XSIZFZ; instead tau=-xsi'*Zg/XSIZFZ
% 980826  hkh  Changed name to qpiSolve
% 980828  hkh  Avoid multiplication if D empty:
%                    if isempty(D), s = []; else s = D(I,:)*p; end
% 980930  hkh  Changed name to qpSolve
% 981001  hkh  Found 3 very serious bugs. Revised the code.
% 981016  hkh  Changed empty Ax to zeros(0,1). Add (:) on input vectors
% 981017  hkh  Nasty bug for test of unconstrained minimum. As W now logical
%              vector test sum(W)==0, not isempty(W), must be done.
% 981018  hkh  Used k on two wrong places (not as counter), changed to kk
%              Infinite loop when adding and deleting same constraint to
%              working set. Added test on number of iterations k
% 981026  hkh  Changed printing levels
% 981028  hkh  Put f_0 = f(x_0) as field in Result
% 981031  hkh  Left and right hand side for equalities got mixed up in phase1
%              calls to lpsimp1 when rhs has different signs
% 981105  hkh  Add code to handle LP, i.e. F empty
% 981107  hkh  Adding comments on Prob.QP field in structure.
%              Sign error in comments on FP problem.
% 981108  hkh  Only use upper bounds from free variables. Use optim2def.
% 981111  hkh  Change FP solution to call lpSimplex
% 981115  hkh  Change alpha == Inf to use isinf.
% 981116  hkh  Improve handling of empty(F)
% 981117  hkh  Skip old transformation method for phase 1.
% 981120  hkh  Return logical basis vector B in QP.B to allow for restart
% 981128  hkh  Change some printing levels for output
% 990827  hkh  Big revision for v2.0 using new lpSimplex features
% 000916  hkh  Adding ExitText
% 001206  hkh  Take real() om eigenvalue/vectors, safeguard to numeric fuzz
% 020822  hkh  Apostrophes missing in save statement
% 040111  hkh  Add call to iniSolve and endSolve
% 040125  hkh  Define field mLin
% 040504  hkh  Avoid FP=1 if unconstrained problem with all variables free
% 040728  med  Pragma removed
% 041126  hkh  Must be some sign error on Lagrange multipliers for linear
%              constraints. Algorithm may accept solution as optimal even if
%              linear constraints are not fulfilled.
%              Made FIX to change sign on lambda if all lambda >=0 and some
%              constraint is not fulfilled
% 041126  hkh  If f is increasing ONLY take backward step IF previous point
%              is feasible.
% 041202  hkh  Added correct computation of Result.xState before every return
% 061028  ango Fix ProbQP.QP.c to work with nonlinear subsolvers
% 070726  frhe Fixed bug introduced when changing all | to |.
% 070907  hkh  SolverFP picked from list, 1st with license, avoid GetSolver
% 080606  med  Switched to iniSolveMini
% 080607  hkh  Use endSolveMini,lp-minos; tomRun, not tomSolve. Comment DEBUG.
% 080607  med  Switched to tomRunMini
% 090717  med  f calculations updated
