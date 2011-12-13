function [ok, s] = releaseverification()
% Run a number of tests and compare results to previously stored values.

% The list of tests is based on releasetest.m, but excluding some tests
% that did not call tomRun, and excluding the quickquide and opt toolbox
% tests.
% Because each call to tomrun gets its own number, numbering will not match
% that of releasetest.

s = struct('run', repmat(struct('probFile',[],'probNo',[],'solName',[], ...
      'preSolve',[],'doSolve',[],'Prob',[],'R',[]),0,0));

s = stdTest(s, 'cls_prob',   20, 'clsSolve');
s = stdTest(s, 'con_prob',   10, 'conSolve');
s = stdTest(s, 'exp_prob',   16, 'expSolve');
mf100 = 'Prob.optParam.MaxFunc = 100';
s = stdTest(s, 'glb_prob',   18, 'glbDirect', mf100);
s = stdTest(s, 'glb_prob',   20, 'glbFast',   mf100);
s = stdTest(s, 'glb_prob',   20, 'glbSolve',  mf100);
s = stdTest(s, 'glc_prob',    2, 'glcCluster', 'Prob.optParam.MaxFunc = 300; Prob.GO.maxFunc1 = 100');
s = stdTest(s, 'glc_prob',    4, 'glcDirect', mf100);
s = stdTest(s, 'glc_prob',    2, 'glcFast',   mf100);
s = stdTest(s, 'glc_prob',    2, 'glcSolve',  mf100);
s = stdTest(s, 'goals_prob',  2, 'goalSolve', mf100); 
s = stdTest(s, 'various_prob',1, 'L1Linsolve'); 
s = stdTest(s, 'lls_prob',    1, 'L1Solve');
s = stdTest(s, 'lp_prob',     5, 'lpSimplex');
s = stdTest(s, 'lls_prob',    1, 'lsei');
s = stdTest(s, 'mip_prob',    6, 'milpSolve');
s = stdTest(s, 'mip_prob',    4, 'mipSolve', 'Prob.optParam.IterPrint = 0');
s = stdTest(s, 'con_prob',    5, 'nlpSolve');
s = stdTest(s, 'qp_prob',     1, 'pdco');
s = stdTest(s, 'qp_prob',     1, 'pdsco', ...
    'Prob.x_U = inf*ones(length(Prob.x_U),1); Prob.x_L = zeros(length(Prob.x_L),1)');
s = stdTest(s, 'qp_prob',     6, 'QLD');
s = stdTest(s, 'qp_prob',     8, 'qpSolve');
s = stdTest(s, 'cls_prob',   11, 'slsSolve');
s = stdTest(s, 'cls_prob',   15, 'sTrustr');
%These don't call tomRun and should probably be in a separate test.
%s = stdTest(s, 'lgo1_prob',   1, 'Tfmin', '', ...
%    '[x, nFunc] = Tfmin(Prob.FUNCS.f, Prob.x_L, Prob.x_U, 1e-6, Prob); R = struct(''x_k'', x, ''FuncEv'', nFunc)');
%s = stdTest(s, 'lgo1_prop', 1, 'Tfzero', 'Prob.FUNCS.f0 = Prob.FUNCS.f', ...
%    ['[xLow, xUpp, ExitFlag] = Tfzero(Prob.x_L, Prob.x_U, Prob, Prob.x_0, 1e-6, 0); ' ...
%     'R = struct(''xLow'', xLow, ''xUpp'', xUpp, ''ExitFlag'', ExitFlag)']);

% Not included because not using tomRun:
% if any(k==Problems)
%     % TLSQR 
%     fprintf('%i: BASE: TESTING TLSQR. \n', k);
%     if(ListOnly==0)
%         N=10^2;
%         A=sparse(N,N);
%         L = N;
%         main = sparse(ones(L,1));
%         off  = sparse(ones(L-1,1));
%         A = diag(main) + diag(off,1) + diag(off,-1);
%         b=rand(N,1);
%         [x] = Tlsqr( N, N, A, [], [], b);
%     end
% end
% k = k + 1;
% 
% if any(k==Problems)
%     % TNNLS
%     fprintf('%i: BASE: TESTING TNLLS. \n', k);
%     if(ListOnly==0)
%         N=10^2;
%         A=sparse(N,N);
%         L = N;
%         main = sparse(ones(L,1));
%         off  = sparse(ones(L-1,1));
%         A = diag(main) + diag(off,1) + diag(off,-1);
%         b=rand(N,1);
%         [x, rNorm, mode, Iter, w] = Tnnls ( full(A), b);
%     end
% end
% k = k + 1;
s = stdTest(s, 'uc_prob',     8, 'ucSolve');
s = stdTest(s, 'cls_prob',   22, 'minos');
s = stdTest(s, 'lp_prob',    10, 'lpopt');
s = stdTest(s, 'qp_prob',     9, 'qpopt');
s = stdTest(s, 'con_prob',    9, 'npsol');
s = stdTest(s, 'cls_prob',    9, 'nlssol');
s = stdTest(s, 'lls_prob',    1, 'lssol');
s = stdTest(s, 'qp_prob',    12, 'sqopt');
s = stdTest(s, 'chs_prob',   85, 'snopt');
  
% s = stdTest(s, '', , '');
 
ok = 1; 
 
% if any(k==Problems)
%   % SNOPT
%   fprintf('%i: SOL: TESTING SNOPT. \n', k);
%   if(ListOnly==0)
%     Prob = probInit('chs_prob', 85);
%     R = tomRun('snopt', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems)
%   % RBFSOLVE
%   fprintf('%i: CGO: TESTING RBFSOLVE. \n', k);
%   if(ListOnly==0)
%     Prob = probInit('glc_prob', 1);
%     Prob.optParam.MaxIter = 5;
%     Prob.optParam.MaxFunc = 10;
%     Prob.optParam.IterPrint = 0;
%     Prob.PriLevOpt = 0;
%     R = tomRun('rbfSolve', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems)
%   % EGO
%   fprintf('%i: CGO: TESTING EGO. \n', k);
%   if(ListOnly==0)
%     Prob = probInit('glc_prob', 1);
%     Prob.optParam.MaxIter = 5;
%     Prob.optParam.MaxFunc = 10;
%     Prob.optParam.IterPrint = 0;
%     Prob.PriLevOpt = 0;
%     R = tomRun('ego', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% % if any(k==Problems) & (plpc | plglx | plsun)
% %   % Xpress
% %   fprintf('%i: Xpress: TESTING Xpress. \n', k);
% %   if(ListOnly==0)
% %     fprintf('%i: Xpress: TESTING LP. \n', k);
% %     Prob = probInit('lp_prob', 6);
% %     R = tomRun('xpress', Prob, PriLev);
% %     fprintf('%i: Xpress: TESTING MILP. \n', k);
% %     Prob = probInit('mip_prob', 4);
% %     R = tomRun('xpress', Prob, PriLev);
% %     fprintf('%i: Xpress: TESTING QP. \n', k);
% %     Prob = probInit('qp_prob', 9);
% %     R = tomRun('xpress', Prob, PriLev);
% %     fprintf('%i: Xpress: TESTING MIQP. \n', k);
% %     Prob = probInit('miqp_prob', 1);
% %     R = tomRun('xpress', Prob, PriLev);
% %     f_k = [f_k; R.f_k];
% %   end
% % end
% % k = k + 1;
% 
% if any(k==Problems) & ~plosx
%   % CPLEX
%   fprintf('%i: CPLEX: TESTING CPLEX. \n', k);
%   if(ListOnly==0)
%     fprintf('%i: CPLEX: TESTING LP. \n', k);
%     Prob = probInit('lp_prob', 6);
%     R = tomRun('cplex', Prob, PriLev);
%     fprintf('%i: CPLEX: TESTING MILP. \n', k);
%     Prob = probInit('mip_prob', 4);
%     R = tomRun('cplex', Prob, PriLev);
%     fprintf('%i: CPLEX: TESTING QP. \n', k);
%     Prob = probInit('qp_prob', 9);
%     R = tomRun('cplex', Prob, PriLev);
%     fprintf('%i: CPLEX: TESTING MIQP. \n', k);
%     Prob = probInit('miqp_prob', 1);
%     R = tomRun('cplex', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems) 
%   % BQPD
%   fprintf('%i: MINLP: TESTING BQPD. \n', k);
%   if(ListOnly==0)
%     fprintf('%i: MINLP: TESTING DENSE BQPD. \n', k);
%     Prob = probInit('qp_prob', 8);
%     R = tomRun('bqpd', Prob, PriLev);
%     Prob.LargeScale = 1;
%     fprintf('%i: MINLP: TESTING SPARSE BQPD. \n', k);
%     R = tomRun('bqpd', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems) 
%   % MIQPBB
%   fprintf('%i: MINLP: TESTING MIQPBB. \n', k);
%   if(ListOnly==0)
%     fprintf('%i: MINLP: TESTING DENSE MIQPBB. \n', k);
%     Prob = probInit('miqp_prob', 2);
%     R = tomRun('miqpBB', Prob, PriLev);
%     Prob.LargeScale = 1;
%     fprintf('%i: MINLP: TESTING SPARSE MIQPBB. \n', k);
%     R = tomRun('miqpBB', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems) 
%   % FILTERSQP
%   fprintf('%i: MINLP: TESTING FILTERSQP. \n', k);
%   if(ListOnly==0)
%     fprintf('%i: MINLP: TESTING DENSE FILTERSQP. \n', k);
%     Prob = probInit('chs_prob', 32);
%     R = tomRun('filterSQP', Prob, PriLev);
%     Prob.LargeScale = 1;
%     Prob.ConsPattern = estConsPattern(Prob);
%     Prob.HessPattern = estHessPattern(Prob);
%     fprintf('%i: MINLP: TESTING SPARSE FILTERSQP. \n', k);
%     R = tomRun('filterSQP', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems) 
%   % MINLPBB
%   fprintf('%i: MINLP: TESTING MINLPBB. \n', k);
%   if(ListOnly==0)
%     fprintf('%i: MINLP: TESTING DENSE MINLPBB. \n', k);
%     Prob = probInit('minlp_prob', 13);
%     R = tomRun('minlpBB', Prob, PriLev);
%     Prob.LargeScale = 1;
%     Prob.HessPattern = estHessPattern(Prob);
%     fprintf('%i: MINLP: TESTING SPARSE MINLPBB. \n', k);
%     R = tomRun('minlpBB', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems) & (plpc | plglx | plsun)
%   % PENSDP
%   fprintf('%i: PENSDP: TESTING PENSDP. \n', k);
%   if(ListOnly==0)
%     Prob = probInit('sdp_prob', 1);
%     R = tomRun('pensdp', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems) & (plpc | plglx | plsun)
%   % PENBMI
%   fprintf('%i: PENBMI: TESTING PENBMI. \n', k);
%   if(ListOnly==0)
%     Prob = probInit('bmi_prob', 1);
%     R = tomRun('penbmi', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems)
%   % KNITRO
%   fprintf('%i: KNITRO: TESTING KNITRO. \n', k);
%   if(ListOnly==0)
%     fprintf('%i: KNITRO: TESTING ALG = 1. \n', k);
%     Prob = probInit('con_prob', 13);
%     Prob.KNITRO.options.ALG = 1;
%     R = tomRun('knitro', Prob, PriLev);
%     Prob.KNITRO.options.ALG = 2;
%     fprintf('%i: KNITRO: TESTING ALG = 2. \n', k);
%     R = tomRun('knitro', Prob, PriLev);
%     Prob.KNITRO.options.ALG = 3;
%     fprintf('%i: KNITRO: TESTING ALG = 3. \n', k);
%     R = tomRun('knitro', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems) & (plpc | plglx)
%   % OQNLP
%   fprintf('%i: OQNLP: TESTING OQNLP. \n', k);
%   if(ListOnly==0)
%     Prob = probInit('glc_prob', 3);
%     Prob.Warning=0;
%     Prob.OQNLP.options.STAGE1_ITERATIONS = 10;
%     Prob.OQNLP.options.MAXTIME = 2;
%     R = tomRun('oqnlp', Prob, PriLev);
%     R = tomRun('msnlp', Prob, PriLev);
%     R = tomRun('lsgrg2', Prob, PriLev);
%     Prob.LargeScale = 1;
%     fprintf('%i: OQNLP: TESTING SPARSE OQNLP. \n', k);
%     R = tomRun('oqnlp', Prob, PriLev);
%     fprintf('%i: OQNLP: TESTING SPARSE MSNLP. \n', k);
%     R = tomRun('msnlp', Prob, PriLev);
%     fprintf('%i: OQNLP: TESTING SPARSE LSGRG2. \n', k);
%     R = tomRun('lsgrg2', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems)
%     % LSGRG2 and MSNLP
%     fprintf('%i: OQNLP: TESTING OQNLP. \n', k);
%     if(ListOnly==0)
%         Prob = probInit('glc_prob', 3);
%         Prob.Warning=0;
%         Prob.OQNLP.options.STAGE1_ITERATIONS = 10;
%         Prob.OQNLP.options.MAXTIME = 2;
%         R = tomRun('msnlp', Prob, PriLev);
%         R = tomRun('lsgrg2', Prob, PriLev);
%         Prob.LargeScale = 1;
%         fprintf('%i: OQNLP: TESTING SPARSE MSNLP. \n', k);
%         R = tomRun('msnlp', Prob, PriLev);
%         fprintf('%i: OQNLP: TESTING SPARSE LSGRG2. \n', k);
%         R = tomRun('lsgrg2', Prob, PriLev);
%         f_k = [f_k; R.f_k];
%     end
% end
% k = k + 1;
% 
% if any(k==Problems) & (plpc | plglx | plsun )
%   % CONOPT
%   fprintf('%i: CONOPT: TESTING CONOPT. \n', k);
%   if(ListOnly==0)
%     Prob = probInit('con_prob', 3);
%     Prob.Warning=0;
%     R = tomRun('conopt', Prob, PriLev);
%     fprintf('%i: CONOPT: TESTING CONOPT - 2. \n', k);
%     Prob.LargeScale = 1;
%     R = tomRun('conopt', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems) & (plpc)
%   % XA
%   fprintf('%i: XA: TESTING XA. \n', k);
%   if(ListOnly==0)
%     fprintf('%i: XA: TESTING LP. \n', k);
%     Prob = probInit('lp_prob', 3);
%     R = tomRun('xa', Prob, PriLev);
%     fprintf('%i: XA: TESTING MILP. \n', k);
%     Prob = probInit('mip_prob', 5);
%     R = tomRun('xa', Prob, PriLev);
%     fprintf('%i: XA: TESTING QP. \n', k);
%     Prob = probInit('qp_prob', 1);
%     R = tomRun('xa', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems) 
%   % NLPQL
%   fprintf('%i: NLPQL: TESTING NLPQLP, NLPJOB, DFNLP. \n', k);
%   if(ListOnly==0)
%     fprintf('%i: NLPQL: TESTING NLPQLP. \n', k);
%     Prob = probInit('con_prob', 3);
%     R = tomRun('nlpqlp', Prob, PriLev);
%     fprintf('%i: NLPQL: TESTING NLPJOB. \n', k);
%     Prob = probInit('mco_prob', 3);
%     R = tomRun('nlpjob', Prob, PriLev);
%     fprintf('%i: NLPQL: TESTING DFNLP. \n', k);
%     Prob = probInit('goals_prob', 1);
%     R = tomRun('dfnlp', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems) 
%   % LGO
%   fprintf('%i: LGO: TESTING LGO. \n', k);
%   if(ListOnly==0)
%     Prob = probInit('lgo2_prob', 22);
%     fprintf('%i: LGO: TESTING OPMODE 0. \n', k);
%     Prob.LGO.options.timelimit = 3;
%     Prob.LGO.options.opmode = 0;
%     R = tomRun('lgo', Prob, PriLev);
%     fprintf('%i: LGO: TESTING OPMODE 1. \n', k);
%     Prob.LGO.options.opmode = 1;
%     R = tomRun('lgo', Prob, PriLev);
%     fprintf('%i: LGO: TESTING OPMODE 2. \n', k);
%     Prob.LGO.options.opmode = 2;
%     R = tomRun('lgo', Prob, PriLev);
%     fprintf('%i: LGO: TESTING OPMODE 3. \n', k);
%     Prob.LGO.options.opmode = 3;
%     R = tomRun('lgo', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems)
%   % GP
%   fprintf('%i: GP: TESTING COPL_GP. \n', k);
%   if(ListOnly==0)
%     Prob = probInit('gp_prob', 12);
%     R = tomRun('gp', Prob, PriLev);
%     f_k = [f_k; R.f_k];
%   end
% end
% k = k + 1;
% 
% if any(k==Problems)
%   % AMPL
%   fprintf('%i: AMPL: TESTING AMPL. \n', k);
%   if(ListOnly==0)
%     a = cd;
%     tomlabdir = which('tomlablic');
%     TOM=fileparts(tomlabdir);
%     cd([TOM '/testprob/nl']);
%     ampltest;
%     cd(a);
%   end
% end
% k = k + 1;
% 
% if any(k==Problems)
%    % GENO
%    fprintf('%i: GENO: TESTING GENO. \n', k);
%    if(ListOnly==0)
%       Prob = probInit('glb_prob', 1);
%       Prob.optParam.MaxIter = 30;
%       R = tomRun('geno', Prob, PriLev);
%       f_k = [f_k; R.f_k];
%    end
% end
% k = k + 1;
% 
% % if any(k==Problems)
% %   DIDO
% %   fprintf('%i: DIDO: TESTING DIDO. \n', k);
% %   if(ListOnly==0)
% %     a = cd;
% %     save a a k;
% %     tomlabdir = which('tomlablic');
% %     TOM=fileparts(tomlabdir);
% %     cd([TOM '/dido']);
% %     startup;
% %     cd('ExampleProblems/MixedProblem');
% %     disp('This could take a couple of minutes...');
% %     J1Main;
% %     load a;
% %     cd(a);
% %   end
% % end
% % k = k + 1;
% 
% % QUICKGUIDE FILES
% 
% if any(k==Problems)
%   fprintf('%i: QUICKGUIDE - RUNNING ALL EXAMPLES. \n', k);
%   if(ListOnly==0)
%     lpQG;
%     mipQG;
%     qpQG;
%     miqpQG;
%     miqqQG;
%     nlpQG;
%     lpconQG;
%     qpconQG;
%     minlpQG;
%     llsQG;
%     millsQG;
%     nllsQG;
%     glbQG;
%     glcQG;
%     sdpQG;
%     bmiQG;
%     minimaxQG;
%     minimaxlinQG;
%     L1QG;
%     L1LinQG;
%     linratQG;
%     goalsQG;
%     simQG;
%     gpQG;
%     expQG;
%     lcpQG;
%     qcpQG;
%     mcpQG;
%     v = version;
%     if str2num(v(1:3)) > 6.49
%         madQG;
%     end
%   end
% end
% k = k +1;
% 
% if any(k==Problems)
%   fprintf('%i: OPT TLBX INTERFACES - RUNNING ALL TESTS. \n', k);
%   if(ListOnly==0)
%     testbintprog;
%     testfmincon;
%     testfminunc;
%     testfgoalattain;
%     testlinprog;
%     testlsqcurvefit;
%     testlsqlin;
%     testlsqnonlin;
%     testlsqnonneg;
%     testquadprog;
%   end
% end
% k = k +1;

% ---- Standard test function
function s = stdTest(s, probFile, probNo, solName, varargin)
  if(length(varargin)>0) && ~isempty(varargin{1})
      preSolve = [varargin{1} ';'];
  else
      preSolve = '';
  end
  if(length(varargin)>1)
      doSolve = [varargin{2} ';'];
  else
      doSolve = '';
  end
  
  fprintf('%i: Testing %s on %s %d. \n', length(s.run)+1, upper(solName), probFile, probNo);
  Prob = probInit(probFile, probNo);
  eval(preSolve);
  if(isempty(doSolve));
    R = tomRun(solName, Prob, -1);
    R = tomRun(solName, Prob, -1); % Repeat with cached functions, as exectution time is logged in R.
  else
    eval(doSolve); % User is responsible for supplying code that sets R.
  end 
  s.run(end+1) = struct('probFile',probFile,'probNo',probNo,'solName',solName, ...
      'preSolve',preSolve,'doSolve',doSolve,'Prob',Prob,'R',R);

