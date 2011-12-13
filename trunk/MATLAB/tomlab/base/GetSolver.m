% GetSolver.m
%
% function Solver = GetSolver(Type, LargeScale, Convex)
%
% GetSolver returns the TOMLAB default solver for different problems
% All standard Types of optimization and also type: FP and DLP, see below
%
% The general minimization problem of QP type is:
%
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
% Type:      String with type of problem:
%  'qp'      Quadratic programming, F is nonempty (QP)
%  'lp'      Linear programming, F is empty, c is nonempty, and c~=0 (LP)
%  'fp'      Feasible point (phase 1) linear programming, F and c empty or 0
%  'dlp'     Dual linear programming. A standard LP problem with a
%            dual feasible initial point available.
% Also Type can be any of the following standard types
%  'uc'      Unconstrained optimization
%  'con'     Nonlinear programming (constrained optimization) (NLP)
%  'ls'      Nonlinear least squares (NLLS)
%  'lls'     Linear least squares (LLS)
%  'cls'     Constrained nonlinear least squares
%  'mip'     Mixed-integer (linear) programming (MIP or MILP)
%  'glb'     Global optimization (box-bounded)
%  'glc'     Global optimization (box-bounded, integer and constrained)
%  'miqp'    Mixed-integer quadratic programming (MIQP)
%  'minlp'   Mixed-integer nonlinear programming (MINLP)
%  'sdp'     Semidefinite programming (SDP) - Linear SDP with LMI constraints
%  'bmi'     Linear SDP with BMI constraints (BMI)
%  'exp'     Parameter estimation in exponential models
%  'ode'     Parameter estimation in ODEs
%  'miqq'    Mixed-Integer Quadratic Programming with Quadratic constraints
%  'gp'      Geometric Programming Problems
%  'mco'     Multi-Criteria Optimization
%  'oc'      Optimal control
%  'lcp'     Standard Linear Complementarity Problem (LCP)
%  'mcp'     Polyhedrally constrained variational inequality Problem or
%            Mixed Complementarity Problem(MCP)
%
% LargeScale If the flag Prob.LargeScale > 0 is set, LargeScale is set true
%
% Convex     If Convex > 0, problem is convex, i.e. F is positive semidefinite
%            for QP problems. Only used for type QP.
%
% OUTPUT PARAMETERS
% Solver     String with name of default solver

%  'nts'     Nonlinear Time Series

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.1.0$
% Written Nov 6, 2000.    Last modified Feb 28, 2009.

function Solver = GetSolver(Type, LargeScale, Convex)

if nargin < 3
   Convex = 1;
   if nargin < 2
      LargeScale = 0;
   end
end

[TomV,os,TV] = tomlabVersion;

switch lower(Type)
  case 'qp'
     % Default TOMLAB QP solver
     % QPOPT
     if Convex & LargeScale
        if TV(9)
           Solver = 'CPLEX';
        elseif TV(4)
           Solver = 'sqopt';
        elseif TV(2)
           Solver = 'qp-minos';
        elseif TV(7)
           Solver = 'bqpd';
        elseif TV(11)
           Solver = 'knitro';
%         elseif TV(8)
%            Solver = 'xpress-mp';
        elseif TV(3)
           Solver = 'qpopt';
        elseif TV(12)
           Solver = 'conopt';
        elseif TV(15)
           Solver = 'xa';
        elseif TV(14)
           Solver = 'oqnlp';
        elseif TV(22)
           Solver = 'msnlp';
        elseif TV(1)
           Solver = 'qld';
        elseif TV(32)
           Solver = 'gurobi';
        else
           error('License for TOMLAB /BASE is missing!!!')
        end
     elseif Convex & ~LargeScale
        if TV(9)
           Solver = 'CPLEX';
        elseif TV(3)
           Solver = 'lssol';
        elseif TV(2)
           Solver = 'qpopt';
        elseif TV(7)
           Solver = 'bqpd';
        elseif TV(4)
           Solver = 'sqopt';
        elseif TV(11)
           Solver = 'knitro';
%         elseif TV(8)
%            Solver = 'xpress-mp';
        elseif TV(12)
           Solver = 'conopt';
        elseif TV(15)
           Solver = 'xa';
        elseif TV(14)
           Solver = 'oqnlp';
        elseif TV(22)
           Solver = 'msnlp';
        elseif TV(1)
           Solver = 'qld';
        elseif TV(32)
           Solver = 'gurobi';
        else
           error('License for TOMLAB /BASE is missing!!!')
        end
     elseif ~Convex & LargeScale
        if TV(4)
           Solver = 'snopt';
        elseif TV(2)
           Solver = 'qp-minos';
        elseif TV(7)
           Solver = 'bqpd';
        elseif TV(11)
           Solver = 'knitro';
        elseif TV(12)
           Solver = 'conopt';
        elseif TV(14)
           Solver = 'oqnlp';
        elseif TV(22)
           Solver = 'msnlp';
        elseif TV(15)
           Solver = 'xa';
        elseif TV(3)
           Solver = 'lssol';
        elseif TV(1)
           Solver = 'qpSolve';
        else
           error('License for TOMLAB /BASE is missing!!!')
        end
     else
        if TV(2)
           Solver = 'qpopt';
        elseif TV(4)
           Solver = 'snopt';
        elseif TV(7)
           Solver = 'bqpd';
        elseif TV(11)
           Solver = 'knitro';
        elseif TV(12)
           Solver = 'conopt';
        elseif TV(14)
           Solver = 'oqnlp';
        elseif TV(22)
           Solver = 'msnlp';
        %elseif TV(15)
        %   Solver = 'xa';
        elseif TV(1)
           Solver = 'qpSolve';
        else
           error('License for TOMLAB /BASE is missing!!!')
        end
     end
  case 'lp'
    if TV(9)
       Solver = 'CPLEX';
    elseif TV(2)
       if LargeScale
          %Solver = 'LP-MINOS';
          Solver = 'MINOS';
       else
          Solver = 'MINOS';
          % Solver = 'lpopt'; % Still MINOS is more reliable
       end
    elseif TV(12)
       Solver = 'conopt';
    elseif TV(15)
       Solver = 'xa';
    elseif TV(11)
       Solver = 'knitro';
    elseif TV(14)
       Solver = 'oqnlp';
    elseif TV(22)
       Solver = 'msnlp';
    elseif TV(1)
       Solver = 'qld';
    elseif TV(32)
       Solver = 'gurobi';
    else
       error('License for TOMLAB /BASE is missing!!!')
    end
  case 'fp'
    if TV(9)
       Solver = 'CPLEX';
    elseif TV(2)
       %Solver = 'LP-MINOS';
       Solver = 'MINOS';
%     elseif TV(8)
%        Solver = 'xpress-mp';
    elseif TV(12)
       Solver = 'conopt';
    elseif TV(15)
       Solver = 'xa';
    elseif TV(11)
       Solver = 'knitro';
    elseif TV(14)
       Solver = 'oqnlp';
    elseif TV(22)
       Solver = 'msnlp';
    elseif TV(1)
       Solver = 'qld';
    else
       error('License for TOMLAB /BASE is missing!!!')
    end
  case 'dlp'
    if TV(2)
       if LargeScale
          %Solver = 'LP-MINOS';
          Solver = 'MINOS';
       else
          Solver = 'MINOS';
          % Solver = 'lpopt'; % Still MINOS is more reliable
       end
    elseif TV(9)
       Solver = 'CPLEX';
%     elseif TV(8)
%        Solver = 'xpress-mp';
    elseif TV(12)
       Solver = 'conopt';
    elseif TV(15)
       Solver = 'xa';
    elseif TV(11)
       Solver = 'knitro';
    elseif TV(14)
       Solver = 'oqnlp';
    elseif TV(22)
       Solver = 'msnlp';
    elseif TV(1)
       Solver = 'qld';
    else
       error('License for TOMLAB /BASE is missing!!!')
    end
  case {'con','lpcon','qpcon'}
    if LargeScale
        if TV(4)
           Solver = 'snopt';
        elseif TV(11)
           Solver = 'knitro';
        elseif TV(12)
           Solver = 'conopt';
        elseif TV(7)
           Solver = 'filterSQP';
        elseif TV(14)
           Solver = 'oqnlp';
        elseif TV(22)
           Solver = 'msnlp';
        elseif TV(3)
           Solver = 'minos';
        elseif TV(2)
           Solver = 'minos';
        elseif TV(1)
           Solver = 'conSolve';
        else
           error('License for TOMLAB /BASE is missing!!!')
        end
    else
        if TV(3)
           Solver = 'npsol';
        elseif TV(4)
           Solver = 'snopt';
        elseif TV(11)
           Solver = 'knitro';
        elseif TV(12)
           Solver = 'conopt';
        elseif TV(7)
           Solver = 'filterSQP';
        elseif TV(14)
           Solver = 'oqnlp';
        elseif TV(22)
           Solver = 'msnlp';
        elseif TV(3)
           Solver = 'minos';
        elseif TV(2)
           Solver = 'minos';
        elseif TV(1)
           Solver = 'conSolve';
        else
           error('License for TOMLAB /BASE is missing!!!')
        end
    end
  case 'uc'
    if LargeScale
        if TV(11)
           Solver = 'knitro';
        elseif TV(4)
           Solver = 'snopt';
        elseif TV(12)
           Solver = 'conopt';
        elseif TV(14)
           Solver = 'oqnlp';
        elseif TV(22)
           Solver = 'msnlp';
        elseif TV(3)
           Solver = 'minos';
        elseif TV(2)
           Solver = 'minos';
        elseif TV(7)
           Solver = 'filterSQP';
        elseif TV(1)
           Solver = 'ucSolve';
        else
           error('License for TOMLAB /BASE is missing!!!')
        end
    else
        if TV(11)
           Solver = 'knitro';
        elseif TV(3)
           Solver = 'npsol';
        elseif TV(4)
           Solver = 'snopt';
        elseif TV(12)
           Solver = 'conopt';
        elseif TV(14)
           Solver = 'oqnlp';
        elseif TV(22)
           Solver = 'msnlp';
        elseif TV(2)
           Solver = 'minos';
        elseif TV(7)
           Solver = 'filterSQP';
        elseif TV(1)
           Solver = 'ucSolve';
        else
           error('License for TOMLAB /BASE is missing!!!')
        end
    end

  case {'ls','cls','exp'}
    if LargeScale
        if TV(4)
           Solver = 'slsSolve';
        elseif TV(3)
           Solver = 'minos';
        elseif TV(2)
           Solver = 'slsSolve';
        elseif TV(11)
           Solver = 'knitro';
        elseif TV(12)
           Solver = 'conopt';
        elseif TV(14)
           Solver = 'oqnlp';
        elseif TV(22)
           Solver = 'msnlp';
        elseif TV(2)
           Solver = 'minos';
        elseif TV(7)
           Solver = 'filterSQP';
        elseif TV(1)
           Solver = 'clsSolve';
        else
           error('License for TOMLAB /BASE is missing!!!')
        end
    else
        if TV(3)
           Solver = 'nlssol';
        elseif TV(4)
           Solver = 'snopt';
        elseif TV(11)
           Solver = 'knitro';
        elseif TV(12)
           Solver = 'conopt';
        elseif TV(14)
           Solver = 'oqnlp';
        elseif TV(22)
           Solver = 'msnlp';
        elseif TV(2)
           Solver = 'minos';
        elseif TV(7)
           Solver = 'filterSQP';
        elseif TV(1)
           Solver = 'clsSolve';
        else
           error('License for TOMLAB /BASE is missing!!!')
        end
    end
  case {'lls'}
    if LargeScale
%         if TV(4)
%            Solver = 'sqopt';
        if TV(4)
           Solver = 'slsSolve';
        elseif TV(11)
           Solver = 'knitro';
        elseif TV(12)
           Solver = 'conopt';
        elseif TV(14)
           Solver = 'oqnlp';
        elseif TV(22)
           Solver = 'msnlp';
        elseif TV(7)
           Solver = 'filterSQP';
%         elseif TV(15)
%            Solver = 'xa';
        elseif TV(3)
           Solver = 'clsSolve'; % clsSolve instead of MINOS
           % Solver = 'minos';
        elseif TV(2)
           Solver = 'clsSolve'; % clsSolve instead of MINOS
           % Solver = 'minos';
        elseif TV(1)
           Solver = 'clsSolve';
        else
           error('License for TOMLAB /BASE is missing!!!')
        end
    else
        if TV(3)
           Solver = 'lssol';
        elseif TV(1)
           Solver = 'lsei'; % Will be picked if not LSSOL present
        elseif TV(11)
           Solver = 'knitro';
        elseif TV(12)
           Solver = 'conopt';
        elseif TV(14)
           Solver = 'oqnlp';
        elseif TV(22)
           Solver = 'msnlp';
        elseif TV(7)
           Solver = 'bqpd';
        elseif TV(15)
           Solver = 'xa';
        elseif TV(4)
           Solver = 'sqopt';
        elseif TV(2)
           Solver = 'minos';
        elseif TV(1)
           Solver = 'clsSolve';
        else
           error('License for TOMLAB /BASE is missing!!!')
        end
    end
  case {'mip'}
    if TV(9)
       Solver = 'CPLEX';
    elseif TV(15)
       Solver = 'xa';
    elseif TV(1)
       Solver = 'mipSolve';
    elseif TV(32)
       Solver = 'gurobi';
    else
       error('License for TOMLAB /BASE is missing!!!')
    end
  case {'glb'}
    if TV(1)
       Solver = 'glbDirect'; %switched from glbFast
    else
       error('License for TOMLAB /BASE is missing!!!')
    end
  case {'glc'}
    if TV(1)
       Solver = 'glcCluster';
    else
       error('License for TOMLAB /BASE is missing!!!')
    end
  case {'miqp'}
    if TV(9)
       Solver = 'CPLEX';
       % elseif TV(8)
       %   Solver = 'xpress-mp';
    elseif TV(7)
       Solver = 'miqpbb';
    elseif TV(32)
       Solver = 'gurobi';
    else
       error('License for TOMLAB /MINLP is missing!!!')
    end
  case {'miqq'}
    if TV(9)
       Solver = 'CPLEX';
    elseif TV(7)
       Solver = 'miqpbb';
    elseif TV(1)
       Solver = 'glcCluster';
    else
       error('License for TOMLAB /BASE is missing!!!')
    end
  case {'minlp'}
    if TV(7)
       Solver = 'minlpbb';
    elseif TV(1)
       Solver = 'glcCluster';
    else
       error('License for TOMLAB /BASE is missing!!!')
    end
  case {'gp'}
    if TV(25)
       Solver = 'coplgp';
    elseif ~TV(1)
       error('License for TOMLAB /BASE is missing!!!')
    else
       Solver = GetSolver('con', LargeScale, 1);
    end
  case {'mco'}
    if TV(16)
       Solver = 'nlpjob';
    elseif TV(1)
       Solver = 'goalSolve';
    else
       error('License for TOMLAB /BASE is missing!!!')
    end
  case {'lcp'}
    if TV(11) %TV(24)
       Solver = 'knitro'; %From path
    elseif ~TV(1)
       error('License for TOMLAB /BASE is missing!!!')
    else
       Solver = GetSolver('con', LargeScale, 0);
    end
  case {'mcp'}
    if TV(11) %TV(24)
       Solver = 'knitro'; %From path
    elseif ~TV(1)
       error('License for TOMLAB /BASE is missing!!!')
    else
       Solver = GetSolver('con', LargeScale, 0);
    end
  case {'oc'}
       Solver = 'propt';
  case {'sdp'}
    if ~TV(1)
       error('License for TOMLAB /BASE is missing!!!')
    else
       Solver = 'PENSDP';
    end
  case {'bmi'}
    if ~TV(1)
       error('License for TOMLAB /BASE is missing!!!')
    else
       Solver = 'PENBMI';
    end
  case {'ode'}
    if ~TV(1)
       error('License for TOMLAB /BASE is missing!!!')
    else
       Solver = 'modfit';
    end
  otherwise
    disp(Type)
    error('Illegal type of optimization problem')
end

% MODIFICATION LOG
% 010726 hkh Differentiate further between sparse and dense con / uc problems
% 020105 hkh Use glbFast and glcCluster as default
% 020701 hkh Use more general license handling, revise all selections
% 020701 hkh Add miqp, miqq, minlp and sdp types
% 030117 hkh Change type miqq to bmi
% 030211 hkh Change QP solver selection
% 030309 hkh Change LP,DLP,CON. First select SOL, otherwise /MINLP
% 041221 hkh Complete revision, prefer cplex to xpress-mp
% 050602 hkh Add GP and other new types
% 061129 ango Change LP order - prefer CPLEX, Xpress, MINOS acc. to hkh.
% 070906 med Commented xpress-mp since not in demo
% 070906 med Fixed some incorrect defaults, dido, path removed
% 080414 hkh Always check if BASE license is OK
% 081006 rut Add 'lpcon' and 'qpcon' to the list
% 090228 med socs removed
