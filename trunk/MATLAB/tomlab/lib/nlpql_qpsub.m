 function [x, u, inform] = nlpql_qpsub(mEQ, C, d, A, b, xl, xu, x0, PriLev, ...
                              Solver, eps)
% Sizes
m = length(b);
n = size(C,1);
mNE = m - mEQ;

mfix = 2;
nfix = 2;

% Bounds on linear constraints
b_L = -b;
b_U = Inf*ones(m,1);
if mEQ > 0
    b_U(1:mEQ) = -b(1:mEQ);
end

% Make sure x0 is within bounds
x0 = max(xl,x0);
x0 = min(xu,x0);

% Maximum iterations
MaxIter = 40*(m+n);

inform = 0;

switch Solver
   case 1 % qpSolve
      Prob = qpAssign(C, d, A, b_L, b_U, xl, xu, x0, ...
             'NLPQL QP subproblem', [], [], [], [], [], [], []);
      Prob.optParam.epsRank = eps;
      Prob.optParam.MaxIter = MaxIter;
      Prob.optParam.xTol = eps;
      Prob.optParam.bTol = eps;
      Result = qpSolve(Prob);
      x = Result.x_k;
   case 2 % QPOPT
      optPar = -999*ones(30,1);
      optPar(10) = eps;
      optPar(11) = eps;
      optPar(30) = MaxIter;
      [Inform, Iter, iState, Ax, Result.v_k, Obj, x] = ...
      qpopt(C, A, [ xl; b_L ], [xu; b_U ], d, 0, x0, optPar);
   case 3 % qp-minos
      Prob = qpAssign(C, d, A, b_L, b_U, xl, xu, x0, ...
             'NLPQL QP subproblem', [], [], [], [], [], [], []);
      optPar = -999*ones(30,1);
      optPar(10) = eps;
      optPar(11) = eps;
      optPar(30) = MaxIter;
      Prob.SOL.optPar = optPar;
      Prob.optParam.bTol = eps;
      Prob.optParam.xTol = eps;
      Result = minosqpTL(Prob);
      x = Result.x_k;
   case 4 % bqpdd
      Prob = qpAssign(C, d, A, b_L, b_U, xl, xu, x0, ...
             'NLPQL QP subproblem', [], [], [], [], [], [], []);
      optPar = -999*ones(20,1);
      optPar(2) = eps;
      Prob.DUNDEE.optPar = optPar;
      Prob.optParam.bTol = eps;
      Prob.optParam.xTol = eps;
      Result = bqpdTL(Prob);
      x = Result.x_k;
      Result.v_k(find(Result.v_k >= 1e20)) = 0;
   case 5 % sqopt
      Prob = qpAssign(sparse(C), sparse(d), sparse(A), b_L, b_U, xl, xu, [], ...
             'NLPQL QP subproblem', [], [], [], [], [], [], []);
      optPar = -999*ones(63,1);
      optPar(11) = eps;
      optPar(12) = eps;
      optPar(30) = MaxIter;
      Prob.SOL.optPar = optPar;
      Prob.optParam.bTol = eps;
      Prob.optParam.xTol = eps;
      Result = sqoptTL(Prob);
      x = Result.x_k;
      switch Result.Inform
         case 3
            inform = 1;
         case 10,22,40
            inform = 3
         otherwise
            inform = 0;
      end
   case 6 % CPLEX
      Prob = qpAssign(sparse(C), sparse(d), sparse(A), b_L, b_U, xl, xu, [], ...
             'NLPQL QP subproblem', [], [], [], [], [], [], []);
%       optPar = -999*ones(63,1);
%       optPar(11) = eps;
%       optPar(12) = eps;
%       optPar(30) = MaxIter;
%       Prob.SOL.optPar = optPar;
%       Prob.optParam.bTol = eps;
%       Prob.optParam.xTol = eps;
      Result = cplexTL(Prob);
   otherwise
      error('Internal error. Contact support@tomopt.com.');
end

% Solvers return the multipliers in different ways;
% Should be fixed, but for now we have to deal with it
u = zeros(m+n+n,1);
u(1:m) = Result.v_k(n+1:n+m);
switch Solver
   case {2, 3}
      mLix = find(Result.v_k(1:n) > 0);
      mUix = find(Result.v_k(1:n) < 0);
      u(m+mLix) = Result.v_k(mLix);
      u(m+nfix+mUix) = -Result.v_k(mUix);
   case {1, 4}
      mLix = find(Result.v_k(1:n) < 0);
      mUix = find(Result.v_k(1:n) > 0);
      u(m+mLix) = Result.v_k(mLix);
      u(m+nfix+mUix) = Result.v_k(mUix);
   case {5}
      mLix = find(Result.v_k(1:n) > 0);
      mUix = find(Result.v_k(1:n) < 0);
      u(m+mLix) = -Result.v_k(mLix);
      u(m+nfix+mUix) = -Result.v_k(mUix);
   otherwise
      mLix = [];
      mUix = [];
end

% switch Solver
%    case 1
%       solvername = 'qpSolve';
%    case 2
%       solvername = 'QPOPT';
%    case 3
%       solvername = 'qpMinos';
%    case 4
%       solvername = 'bqpdd';
%    case 5
%       solvername = 'sqopt';
%    case 6
%       solvername = 'CPLEX';
% end
% disp(['                                   ' solvername]);
% disp(' ');
