%
% Compute loop 2 in multiMin with parfor using Prob.Threads threads


% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written July 13, 2006.    Last modified July 23, 2011.

function [Eval,C02,CX1,CX2,FX1,EX2,EX3,XX1] = multiMin_2_par(Prob,X0,F0, ...
          IntVars,localSolver,SOLSNOPT,PriLev);
global NLP_x NLP_f
b_L     = Prob.b_L;
b_U     = Prob.b_U;
c_L     = Prob.c_L;
c_U     = Prob.c_U;
dLin    = Prob.mLin;
dNoL    = Prob.mNonLin;
absviol = 0;
[n,M]   = size(X0);
C02     = zeros(1,M);
CX1     = zeros(1,M);
CX2     = zeros(1,M);
FX1     = zeros(1,M);
EX2     = zeros(1,M);
EX3     = zeros(1,M);
XX1     = zeros(n,M);
Ev1     = zeros(M,1);
Ev2     = zeros(M,1);
Ev3     = zeros(M,1);
Ev4     = zeros(M,1);
Ev5     = zeros(M,1);
% Local Optimization for each x_0
parfor P=1:M
    ProbI = Prob;
    x_0   = X0(:,P);
    NLP_x = x_0;
    NLP_f = F0(P);
    if dNoL > 0
       c_0       = nlp_c(x_0, ProbI);
       [NV0 NVi] = consViol(c_0,c_L,c_U,absviol);
       C02(P)   = NV0;
    end
    ProbI.x_0     = x_0;
    ProbI.Name    = [Prob.Name, ' - Trial ' num2str(P)];
    if ~isempty(IntVars)          % Fix integer variables
       ProbI.x_L(IntVars) = x_0(IntVars);
       ProbI.x_U(IntVars) = x_0(IntVars);
    end
    r            = tomRun(localSolver,ProbI,PriLev-7);
    %PrintResult(r,2);
    if SOLSNOPT == 1 & r.FuncEv == 0  % NOT NEEDED: & Inform ~= -999
       Ev1(P)   = r.iwCount(3);
       Ev2(P)   = sum(r.iwCount(4:6));
       Ev3(P)   = r.iwCount(7);
    else
       Ev1(P)   = r.FuncEv;
       Ev2(P)   = r.GradEv;
       Ev3(P)   = r.ConstrEv;
    end
    Ev4(P)      = r.HessEv;
    Ev5(P)      = r.Iter;
    Inform       = r.Inform;
    f_k          = r.f_k;
    FX1(P)       = f_k;
    EX2(P)       = r.ExitFlag;
    EX3(P)       = Inform;
    % Avoid tiny value outside bounds
    x_k          = min(ProbI.x_U,max(ProbI.x_L,r.x_k));
    XX1(:,P)     = x_k;
    if dLin > 0                   % Linear constraint violation
       z         = r.Ax;
       if isempty(z)
          z      = ProbI.A*x_k;
       end
       [LV LVi]  = consViol(z,b_L,b_U,absviol);
       CX1(P)    = LV;
    end
    if dNoL > 0                  % Nonlinear constraint violation
       c_k       = r.c_k;
       if isempty(c_k)
          c_k    = nlp_c(x_k, ProbI);
       end
       [NV NVi]  = consViol(c_k,c_L,c_U,absviol);
       CX2(P)   = NV;
    end
end
Eval(1) = sum(Ev1);
Eval(2) = sum(Ev2);
Eval(3) = sum(Ev3);
Eval(4) = sum(Ev4);
Eval(5) = sum(Ev5);

% MODIFICATION LOG:
%
% 110723  hkh  Written
