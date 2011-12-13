% function releasetest(Platform, PriLev, Problems, ListOnly)
%
% releasetest runs one problem for each solver for TOMLAB Init Files.
%
% INPUT:
% Platform    1: win32
%             2: linux32
%             3: sun,
%             4: macosx 32-bit intel
%             5: linux64 
%             6: win64   
%             7: macosx64 64-bit intel
%
% PriLev      Print level in tomRun routine (0 default).
% Problems    Index of problems to test. Default all.
% ListOnly    List problem indices only. Default 0 (no).
%
% OUTPUT:
% No output as of now

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Nov 4, 2004.  Last modified Sept 10, 2009.

function releasetest(Platform, PriLev, Problems, ListOnly)

global f_k

% A development license is needed to run this file.
if nargin < 4
  ListOnly = 0;
  if nargin < 3
    Problems = (1:100)';
    if nargin < 2
      PriLev=0;
      if nargin < 1
      Platform=[];
      end 
    end 
  end
end

if isempty(Platform)
  switch(computer)
   case 'PCWIN'
    Platform = 1;
    
   case 'PCWIN64'
    Platform = 6;
    
   case 'GLNX86'
    Platform = 2;
    
   case 'SOL2'
    Platform = 3;
    
   case {'MAC', 'MACI'}
    Platform = 4;
   
   case {'MACI64'}
    Platform = 7;
   
   case 'GLNXA64'
    Platform = 5;
    
   otherwise
    error('Unknown platform');
    
  end
end

plpc    = (Platform == 1);
plwin32 = (Platform == 1);
plglx   = (Platform == 2);
plsun   = (Platform == 3);
plosx   = (Platform == 4);
plglx64 = (Platform == 5);
plwin64 = (Platform == 6);
plosx64 = (Platform == 7);

% f_k vector for comparison

f_k = [];

k = 1;

if any(k==Problems)
  % CLSSOLVE
  fprintf('%i: BASE: TESTING CLSSOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('cls_prob', 20);
    R = tomRun('clsSolve', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % CONSOLVE
  fprintf('%i: BASE: TESTING CONSOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('con_prob', 10);
    R = tomRun('conSolve', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
   % EXPSOLVE
   fprintf('%i: BASE: TESTING EXPSOLVE. \n', k);
   if(ListOnly==0)
   Prob = probInit('exp_prob', 16);
   R = tomRun('expSolve', Prob, PriLev);
   f_k = [f_k; R.f_k];
   end
end
k = k + 1;

if any(k==Problems)
  % GLBDIRECT
  fprintf('%i: BASE: TESTING GLBDIRECT. \n', k);
  if(ListOnly==0)
    Prob = probInit('glb_prob', 18);
    Prob.optParam.MaxFunc = 100;
    R = tomRun('glbDirect', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % GLBFAST
  fprintf('%i: BASE: TESTING GLBFAST. \n', k);
  if(ListOnly==0)
    Prob = probInit('glb_prob', 20);
    Prob.optParam.MaxFunc = 100;
    R = tomRun('glbFast', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % GLBSOLVE
  fprintf('%i: BASE: TESTING GLBSOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('glb_prob', 20);
    Prob.optParam.MaxFunc = 100;
    R = tomRun('glbSolve', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % GLCCLUSTER
  fprintf('%i: BASE: TESTING GLCCLUSTER. \n', k);
  if(ListOnly==0)
    Prob = probInit('glc_prob', 2);
    Prob.optParam.MaxFunc = 300;
    Prob.GO.maxFunc1 = 100;
    R = tomRun('glcCluster', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % GLCDIRECT
  fprintf('%i: BASE: TESTING GLCDIRECT. \n', k);
  if(ListOnly==0)
    Prob = probInit('glc_prob', 4);
    Prob.optParam.MaxFunc = 100;
    R = tomRun('glcDirect', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % GLCFAST
  fprintf('%i: BASE: TESTING GLCFAST. \n', k);
  if(ListOnly==0)
    Prob = probInit('glc_prob', 2);
    Prob.optParam.MaxFunc = 100;
    R = tomRun('glcFast', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % GLCSOLVE
  fprintf('%i: BASE: TESTING GLCSOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('glc_prob', 2);
    Prob.optParam.MaxFunc = 100;
    R = tomRun('glcSolve', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems) & 0 
  % GOALSOLVE
  fprintf('%i: BASE: TESTING GOALSOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('goals_prob', 2);
    Prob.optParam.MaxFunc = 100;
    Prob = ProbCheck(Prob, 'goalSolve');
    R = goalSolve(Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
    % L1LINSOLVE
    fprintf('%i: BASE: TESTING L1LINSOLVE. \n', k);
    if(ListOnly==0)
        Name='L1LinSolve test example';
        n = 6;
        x_L = -10*ones(n,1);
        x_U =  10*ones(n,1);
        x_0 = (x_L + x_U) / 2;

        C = spdiags([1 2 3 4 5 6]', 0, n, n);
        y = 1.5*ones(n,1);
        A = [1 1 0 0 0 0]; b_L = 1; b_U = 1;

        Prob = llsAssign(C, y, x_L, x_U, Name, x_0, ...
            [], [], [], ...
            A, b_L, b_U);
 
        Prob.LS.damp = 1;
        Prob.LS.L = spdiags(ones(6,1)*0.01, 0, 6, 6);

        Prob.SolverL1 = 'lpSimplex';
        R = L1LinSolve(Prob, 0);
        f_k = [f_k; R.f_k];
    end
end
k = k + 1;

if any(k==Problems)
  % L1SOLVE
  fprintf('%i: BASE: TESTING L1SOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('lls_prob', 1);
    R = L1Solve(Prob, 0);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % LPSIMPLEX
  fprintf('%i: BASE: TESTING LPSIMPLEX. \n', k);
  if(ListOnly==0)
    Prob = probInit('lp_prob', 5);
    R = tomRun('lpSimplex', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % LSEI
  fprintf('%i: BASE: TESTING LSEI. \n', k);
  if(ListOnly==0)
    Prob = probInit('lls_prob', 1);
    R = tomRun('lsei', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % MILPSOLVE
  fprintf('%i: BASE: TESTING MILPSOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('mip_prob', 6);
    R = tomRun('milpSolve', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % MIPSOLVE
  fprintf('%i: BASE: TESTING MIPSOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('mip_prob', 4);
    Prob.optParam.IterPrint = 0;
    R = tomRun('mipSolve', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % NLPSOLVE
  fprintf('%i: BASE: TESTING NLPSOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('con_prob', 5);
    R = tomRun('nlpSolve', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % PDCO
  fprintf('%i: BASE: TESTING PDCO. \n', k);
  if(ListOnly==0)
    Prob = probInit('qp_prob', 1);
    R = tomRun('pdco', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % PDSCO
  fprintf('%i: BASE: TESTING PDSCO. \n', k);
  if(ListOnly==0)
    Prob = probInit('qp_prob', 1);
    Prob.x_U = inf*ones(length(Prob.x_U),1);
    Prob.x_L = zeros(length(Prob.x_L),1);
    R = tomRun('pdsco', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % QLD
  fprintf('%i: BASE: TESTING QLD. \n', k);
  if(ListOnly==0)
    Prob = probInit('qp_prob', 6);
    R = tomRun('QLD', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % QPSOLVE
  fprintf('%i: BASE: TESTING QPSOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('qp_prob', 8);
    R = tomRun('qpSolve', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % SLSSOLVE
  fprintf('%i: BASE: TESTING SLSSOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('cls_prob', 11);
    Prob = ProbCheck(Prob, 'slsSolve');
    R = slsSolve(Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % STRUSTR
  fprintf('%i: BASE: TESTING STRUSTR. \n', k);
  if(ListOnly==0)
    Prob = probInit('cls_prob', 15);
    R = tomRun('sTrustr', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % TFMIN
  fprintf('%i: BASE: TESTING TFMIN. \n', k);
  if(ListOnly==0)
    Prob = probInit('lgo1_prob', 1);
    [x, nFunc] = Tfmin(Prob.FUNCS.f, Prob.x_L, Prob.x_U, 1e-6, Prob);
    f_k = [f_k; nlp_f(x, Prob)];
  end
end
k = k + 1;

if any(k==Problems)
  % TFZERO
  fprintf('%i: BASE: TESTING TFZERO. \n', k);
  if(ListOnly==0)
    Prob = probInit('lgo1_prob', 1);
    Prob.FUNCS.f0 = Prob.FUNCS.f;
    [xLow, xUpp, ExitFlag] = Tfzero(Prob.x_L, Prob.x_U, Prob, Prob.x_0, 1e-6, 0);
    f_k = [f_k; nlp_f(xLow, Prob)];
  end
end
k = k + 1;

if any(k==Problems)
    % TLSQR
    fprintf('%i: BASE: TESTING TLSQR. \n', k);
    if(ListOnly==0)
        N=10^2;
        A=sparse(N,N);
        L = N;
        main = sparse(ones(L,1));
        off  = sparse(ones(L-1,1));
        A = diag(main) + diag(off,1) + diag(off,-1);
        b=rand(N,1);
        [x] = Tlsqr( N, N, A, [], [], b);
        if PriLev>0
            fprintf('|Ax-b|=%16.8g, |inv(A)*b-x| = %16.8g\n',norm(A*x-b),norm(inv(A)*b-x))
        end
    end
end
k = k + 1;

if any(k==Problems)
    % TNNLS
    fprintf('%i: BASE: TESTING TNLLS. \n', k);
    if(ListOnly==0)
        N=10^2;
        A=sparse(N,N);
        L = N;
        main = sparse(ones(L,1));
        off  = sparse(ones(L-1,1));
        A = diag(main) + diag(off,1) + diag(off,-1);
        b=rand(N,1);
        [x, rNorm, mode, Iter, w] = Tnnls ( full(A), b);
        if PriLev>0
            fprintf('|Ax-b|=%16.8g, min(x)=%16.8g\n',norm(A*x-b),min(x));
        end
    end
end
k = k + 1;

if any(k==Problems)
  % UCSOLVE
  fprintf('%i: BASE: TESTING UCSOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('uc_prob', 8);
    R = tomRun('ucSolve', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % MINOS
  fprintf('%i: SOL: TESTING MINOS. \n', k);
  if(ListOnly==0)
    Prob = probInit('cls_prob', 22);
    R = tomRun('minos', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % LPOPT
  fprintf('%i: SOL: TESTING LPOPT. \n', k);
  if(ListOnly==0)
    Prob = probInit('lp_prob', 10);
    R = tomRun('lpopt', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % QPOPT
  fprintf('%i: SOL: TESTING QPOPT. \n', k);
  if(ListOnly==0)
    Prob = probInit('qp_prob', 9);
    R = tomRun('qpopt', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % NPSOL
  fprintf('%i: SOL: TESTING NPSOL. \n', k);
  if(ListOnly==0)
    Prob = probInit('con_prob', 9);
    R = tomRun('npsol', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % NLSSOL
  fprintf('%i: SOL: TESTING NLSSOL. \n', k);
  if(ListOnly==0)
    Prob = probInit('cls_prob', 9);
    R = tomRun('nlssol', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % LSSOL
  fprintf('%i: SOL: TESTING LSSOL. \n', k);
  if(ListOnly==0)
    Prob = probInit('lls_prob', 1);
    R = tomRun('lssol', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % SQOPT
  fprintf('%i: SOL: TESTING SQOPT. \n', k);
  if(ListOnly==0)
    Prob = probInit('qp_prob', 12);
    R = tomRun('sqopt', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % SNOPT
  fprintf('%i: SOL: TESTING SNOPT. \n', k);
  if(ListOnly==0)
    Prob = probInit('chs_prob', 85);
    R = tomRun('snopt', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % snoptA
  fprintf('%i: SOL: TESTING SNOPTA. \n', k);
  if(ListOnly==0)
    snoptA_Ex1;
  end
end
k = k + 1;

if any(k==Problems)
  % RBFSOLVE
  fprintf('%i: CGO: TESTING RBFSOLVE. \n', k);
  if(ListOnly==0)
    Prob = probInit('glc_prob', 1);
    Prob.optParam.MaxIter = 5;
    Prob.optParam.MaxFunc = 10;
    Prob.optParam.IterPrint = 0;
    Prob.PriLevOpt = 0;
    R = tomRun('rbfSolve', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % EGO
  fprintf('%i: CGO: TESTING EGO. \n', k);
  if(ListOnly==0)
    Prob = probInit('glc_prob', 1);
    Prob.optParam.MaxIter = 5;
    Prob.optParam.MaxFunc = 10;
    Prob.optParam.IterPrint = 0;
    Prob.PriLevOpt = 0;
    R = tomRun('ego', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

% if any(k==Problems) & (plpc | plglx | plsun)
%   % Xpress
%   fprintf('%i: Xpress: TESTING Xpress. \n', k);
%   if(ListOnly==0)
%     fprintf('%i: Xpress: TESTING LP. \n', k);
%     Prob = probInit('lp_prob', 6);
%     R = tomRun('xpress', Prob, PriLev);
%     fprintf('%i: Xpress: TESTING MILP. \n', k);
%     Prob = probInit('mip_prob', 4);
%     R = tomRun('xpress', Prob, PriLev);
%     fprintf('%i: Xpress: TESTING QP. \n', k);
%     Prob = probInit('qp_prob', 9);
%     R = tomRun('xpress', Prob, PriLev);
%     fprintf('%i: Xpress: TESTING MIQP. \n', k);
%     Prob = probInit('miqp_prob', 1);
%     R = tomRun('xpress', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;

if any(k==Problems)
  % CPLEX
  fprintf('%i: CPLEX: TESTING CPLEX. \n', k);
  if(ListOnly==0)
    fprintf('%i: CPLEX: TESTING LP. \n', k);
    Prob = probInit('lp_prob', 6);
    R = tomRun('cplex', Prob, PriLev);
    fprintf('%i: CPLEX: TESTING MILP. \n', k);
    Prob = probInit('mip_prob', 4);
    R = tomRun('cplex', Prob, PriLev);
    fprintf('%i: CPLEX: TESTING QP. \n', k);
    Prob = probInit('qp_prob', 9);
    R = tomRun('cplex', Prob, PriLev);
    fprintf('%i: CPLEX: TESTING MIQP. \n', k);
    Prob = probInit('miqp_prob', 1);
    R = tomRun('cplex', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

try
if any(k==Problems) & ~plosx
  % GUROBI
  fprintf('%i: GUROBI: TESTING GUROBI. \n', k);
  if(ListOnly==0)
    fprintf('%i: GUROBI: TESTING LP. \n', k);
    Prob = probInit('lp_prob', 6);
    R = tomRun('gurobi', Prob, PriLev);
    fprintf('%i: GUROBI: TESTING MILP. \n', k);
    Prob = probInit('mip_prob', 4);
    R = tomRun('gurobi', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
catch
  disp(lasterror)
end

  k = k + 1;


if any(k==Problems) 
  % BQPD
  fprintf('%i: MINLP: TESTING BQPD. \n', k);
  if(ListOnly==0)
    fprintf('%i: MINLP: TESTING DENSE BQPD. \n', k);
    Prob = probInit('qp_prob', 8);
    R = tomRun('bqpd', Prob, PriLev);
    Prob.LargeScale = 1;
    fprintf('%i: MINLP: TESTING SPARSE BQPD. \n', k);
    R = tomRun('bqpd', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems) 
  % MIQPBB
  fprintf('%i: MINLP: TESTING MIQPBB. \n', k);
  if(ListOnly==0)
    fprintf('%i: MINLP: TESTING DENSE MIQPBB. \n', k);
    Prob = probInit('miqp_prob', 2);
    R = tomRun('miqpBB', Prob, PriLev);
    Prob.LargeScale = 1;
    fprintf('%i: MINLP: TESTING SPARSE MIQPBB. \n', k);
    R = tomRun('miqpBB', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems) 
  % FILTERSQP
  fprintf('%i: MINLP: TESTING FILTERSQP. \n', k);
  if(ListOnly==0)
    fprintf('%i: MINLP: TESTING DENSE FILTERSQP. \n', k);
    Prob = probInit('chs_prob', 32);
    R = tomRun('filterSQP', Prob, PriLev);
    Prob.LargeScale = 1;
    Prob.Warning = 0;
    Prob = ProbCheck(Prob,'filterSQP',3,3);
    Prob.ConsPattern = 0;
    Prob = iniSolve(Prob,3,2,1);
    % Testing estConsPattern med 5 trials
    Prob.ConsPattern = estConsPattern(Prob,5);
    % ConsPattern = Prob.ConsPattern
    fprintf('%i: MINLP: TESTING SPARSE FILTERSQP. \n', k);
    R = tomRun('filterSQP', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems) 
  % MINLPBB
  fprintf('%i: MINLP: TESTING MINLPBB. \n', k);
  if(ListOnly==0)
    fprintf('%i: MINLP: TESTING DENSE MINLPBB. \n', k);
    Prob = probInit('minlp_prob', 13);
    R = tomRun('minlpBB', Prob, PriLev);
    Prob.LargeScale = 1;
    Prob = ProbCheck(Prob,'minlpBB',3,3);
    Prob.Warning = 0;
    Prob = iniSolve(Prob,3,2,1);
    fprintf('%i: MINLP: TESTING SPARSE MINLPBB. \n', k);
    R = tomRun('minlpBB', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems) & (plpc | plwin64 | plglx | plsun)
  % PENSDP
  fprintf('%i: PENSDP: TESTING PENSDP. \n', k);
  if(ListOnly==0)
    Prob = probInit('sdp_prob', 1);
    R = tomRun('pensdp', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems) & (plpc | plwin64 | plglx | plsun)
  % PENBMI
  fprintf('%i: PENBMI: TESTING PENBMI. \n', k);
  if(ListOnly==0)
    Prob = probInit('bmi_prob', 1);
    R = tomRun('penbmi', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % KNITRO
  fprintf('%i: KNITRO: TESTING KNITRO. \n', k);
  if(ListOnly==0)
    fprintf('%i: KNITRO: TESTING ALG = 1. \n', k);
    Prob = probInit('con_prob', 13);
    Prob.KNITRO.options.ALG = 1;
    R = tomRun('knitro', Prob, PriLev);
    Prob.KNITRO.options.ALG = 2;
    fprintf('%i: KNITRO: TESTING ALG = 2. \n', k);
    R = tomRun('knitro', Prob, PriLev);
    Prob.KNITRO.options.ALG = 3;
    fprintf('%i: KNITRO: TESTING ALG = 3. \n', k);
    R = tomRun('knitro', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems) & (plwin32 | plwin64 | plglx)
  % OQNLP
  fprintf('%i: OQNLP: TESTING OQNLP. \n', k);
  if(ListOnly==0)
    Prob = probInit('glc_prob', 3);
    Prob.Warning=0;
    Prob.OQNLP.options.STAGE1_ITERATIONS = 10;
    Prob.OQNLP.options.MAXTIME = 2;
    R = tomRun('oqnlp', Prob, PriLev);
    R = tomRun('msnlp', Prob, PriLev);
    R = tomRun('lsgrg2', Prob, PriLev);
    Prob.LargeScale = 1;
    fprintf('%i: OQNLP: TESTING SPARSE OQNLP. \n', k);
    R = tomRun('oqnlp', Prob, PriLev);
    fprintf('%i: OQNLP: TESTING SPARSE MSNLP. \n', k);
    R = tomRun('msnlp', Prob, PriLev);
    fprintf('%i: OQNLP: TESTING SPARSE LSGRG2. \n', k);
    R = tomRun('lsgrg2', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
    % LSGRG2 and MSNLP
    fprintf('%i: OQNLP: TESTING OQNLP. \n', k);
    if(ListOnly==0)
        Prob = probInit('glc_prob', 3);
        Prob.Warning=0;
        Prob.OQNLP.options.STAGE1_ITERATIONS = 10;
        Prob.OQNLP.options.MAXTIME = 2;
        R = tomRun('msnlp', Prob, PriLev);
        R = tomRun('lsgrg2', Prob, PriLev);
        Prob.LargeScale = 1;
        fprintf('%i: OQNLP: TESTING SPARSE MSNLP. \n', k);
        R = tomRun('msnlp', Prob, PriLev);
        fprintf('%i: OQNLP: TESTING SPARSE LSGRG2. \n', k);
        R = tomRun('lsgrg2', Prob, PriLev);
        f_k = [f_k; R.f_k];
    end
end
k = k + 1;

if any(k==Problems) & (plpc | plglx | plsun | plosx64 )
  % CONOPT
  fprintf('%i: CONOPT: TESTING CONOPT. \n', k);
  if(ListOnly==0)
    Prob = probInit('con_prob', 3);
    Prob.Warning=0;
    R = tomRun('conopt', Prob, PriLev);
    fprintf('%i: CONOPT: TESTING CONOPT - 2. \n', k);
    Prob.LargeScale = 1;
    R = tomRun('conopt', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems) 
  % NLPQL
  fprintf('%i: NLPQL: TESTING NLPQLP, NLPJOB, DFNLP. \n', k);
  if(ListOnly==0)
    fprintf('%i: NLPQL: TESTING NLPQLP. \n', k);
    Prob = probInit('con_prob', 3);
    R = tomRun('nlpqlp', Prob, PriLev);
    fprintf('%i: NLPQL: TESTING NLPJOB. \n', k);
    Prob = probInit('mco_prob', 3);
    R = tomRun('nlpjob', Prob, PriLev);
    fprintf('%i: NLPQL: TESTING DFNLP. \n', k);
    Prob = probInit('goals_prob', 1);
    R = tomRun('dfnlp', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems) 
  % LGO
  fprintf('%i: LGO: TESTING LGO. \n', k);
  if(ListOnly==0)
    Prob = probInit('lgo2_prob', 22);
    fprintf('%i: LGO: TESTING OPMODE 0. \n', k);
    Prob.LGO.options.timelimit = 3;
    Prob.LGO.options.opmode = 0;
    R = tomRun('lgo', Prob, PriLev);
    fprintf('%i: LGO: TESTING OPMODE 1. \n', k);
    Prob.LGO.options.opmode = 1;
    R = tomRun('lgo', Prob, PriLev);
    fprintf('%i: LGO: TESTING OPMODE 2. \n', k);
    Prob.LGO.options.opmode = 2;
    R = tomRun('lgo', Prob, PriLev);
    fprintf('%i: LGO: TESTING OPMODE 3. \n', k);
    Prob.LGO.options.opmode = 3;
    R = tomRun('lgo', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % GP
  fprintf('%i: GP: TESTING COPL_GP. \n', k);
  if(ListOnly==0)
    Prob = probInit('gp_prob', 12);
    R = tomRun('gp', Prob, PriLev);
    f_k = [f_k; R.f_k];
  end
end
k = k + 1;

if any(k==Problems)
  % AMPL
  fprintf('%i: AMPL: TESTING AMPL. \n', k);
  if(ListOnly==0)
    a = cd;
    tomlabdir = which('tomlablic');
    TOM=fileparts(tomlabdir);
    cd([TOM '/testprob/nl']);
    ampltest;
    cd(a);
  end
end
k = k + 1;

if any(k==Problems)
   % GENO
   fprintf('%i: GENO: TESTING GENO. \n', k);
   if(ListOnly==0)
      Prob = probInit('glb_prob', 1);
      Prob.optParam.MaxIter = 30;
      R = tomRun('geno', Prob, PriLev);
      f_k = [f_k; R.f_k];
   end
end
k = k + 1;

% QUICKGUIDE FILES

if any(k==Problems)
  fprintf('%i: QUICKGUIDE - RUNNING ALL EXAMPLES. \n', k);
  if(ListOnly==0)
    lpQG;
    mipQG;
    qpQG;
    miqpQG;
    miqqQG;
    nlpQG;
    lpconQG;
    qpconQG;
    minlpQG;
    llsQG;
    millsQG;
    nllsQG;
    glbQG;
    glcQG;
    if ~plosx && ~plosx64
      sdpQG;
      bmiQG;
    end
    minimaxQG;
    minimaxlinQG;
    L1QG;
    L1LinQG;
    linratQG;
    goalsQG;
    simQG;
    gpQG;
    expQG;
    lcpQG;
    qcpQG;
    mcpQG;
    v = version;
    if str2num(v(1:3)) > 6.9
        funchandleQG;
    end
    if ~plosx
      madQG;
    end
    close all
    proptQG;
  end
end
k = k +1;

if any(k==Problems)
  fprintf('%i: OPT TLBX INTERFACES - RUNNING ALL TESTS. \n', k);
  if(ListOnly==0)
    testbintprog;
    testfmincon;
    testfminunc;
    testfgoalattain;
    testlinprog;
    testlsqcurvefit;
    testlsqlin;
    testlsqnonlin;
    testlsqnonneg;
    testquadprog;
  end
end
k = k +1;

% MODIFICATION LOG:
%
% 041104 med  Written.
% 041105 frhe Added parameter: TestList. Removed static problem indices.
%             Added PATH.
% 041123 hkh  Change call to tomRun
% 050112 med  Added MSNLP and LSGRG2
% 050112 med  Quickguide examples added
% 050113 frhe Testlist check added to quickguide examples
% 050113 frhe Enabled LGO and PENOPT solvers and added platform
%             checks for all solvers. TestList => ListOnly
% 050117 med  Bug for LGO removed, prob number was i
% 050117 med  LGO options changed
% 050210 med  lpconQG and qpconQG added
% 050211 med  simQG added
% 050222 med  BARNLP/SPRNLP solvers added
% 050309 med  Changed problem used for BARLS and SPRLS
% 050309 med  Removed PATH
% 050601 med  Added TOMLAB /GP, glbDirect and glcDirect
% 050701 ango Better at finding correct platform
% 060203 med  Added more quickguide problems
% 061201 ango Added GENO
% 070307 ango KNITRO everywhere
% 080604 hkh  Call ProbCheck, iniSolve before estConsPattern, estHessPattern
% 080619 med  PROPT and SOCS now in quickguide folder
% 090228 med  socs removed
% 090910 hkh  New call to Tnnls
