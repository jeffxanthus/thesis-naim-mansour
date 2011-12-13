% nlpSolve.m:
%
% nlpSolve implements the Filter SQP by Roger Fletcher and Sven Leyffer.
%
% Filter SQP: presented in "Nonlinear Programming without a penalty
% function". All page references and notations are taken from this paper.
%
% Minimization problem:
%        min          f(x)
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%              c_L <= c(x) <= c_U
%
% function Result = nlpSolve(Prob, varargin)
%
% nlpSolve runs different methods to obtain the gradient g and the Hessian H,
% dependent on the parameters Prob.Solver.Alg, Prob.NumDiff, Prob.ADObj
% and if a user supplied routine for the Hessian, stored in
% Prob.FUNCS.H is available.
% Solver.Alg=1 gives the Filter SQP
% Solver.Alg=2 gives quasi-Newton BFGS for the Filter SQP
%
% The table gives the different possibilities for Alg=0 and 1
%
% Solver.Alg  NumDiff    ADObj isempty(FUNCS.H)  Hessian computation
%     0          0         0           0         Analytic Hessian
%     0        any       any           1         BFGS
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
% INPUT PARAMETERS
%   Use conAssign.m (or probAssign) to initialize the Prob structure in the
%   Tomlab Quick format.  Fields used in structure Prob:
%
%   x_0      Starting point
%   x_L      Lower bounds for x
%   x_U      Upper bounds for x
%   b_L      Lower bounds for linear constraints
%   b_U      Upper bounds for linear constraints
%   A        Linear constraint matrix
%   c_L      Lower bounds for nonlinear constraints
%   c_U      Upper bounds for nonlinear constraints
%   f_Low    A lower bound on the optimal function value
%            Used in convergence tests, f_k(x_k) <= f_Low. Only a feasible
%            point x_k is accepted
%
%   FUNCS.f  The routine to compute the function as a string
%   FUNCS.g  The routine to compute the gradient as a string
%   FUNCS.H  The routine to compute the Hessian as a string
%   FUNCS.c  The routine to evaluate the constraints as a string
%   FUNCS.dc The routine to compute the gradient of the constraints
%   FUNCS.d2c The routine to compute the second derivatives of the constraints
%   NumDiff  How to obtain derivatives (gradient, Hessian)
%   ConsDiff Differentiation method for the constraint Jacobian
%            0 = analytic, 1-5 different numerical methods
%   SolverQP Name of the solver used for QP subproblems. If empty,
%            picked from a list, best available with a license
%   SolverFP Name of the solver used for FP (feasible point) subproblems.
%            If empty, picked from a list, best available with a license
%
%   PriLevOpt Print level: 0 Silent, 1 Final result, 2 Each iteration, short
%   	     3 Each iteration, more info, 4 Matrix update information
%   optParam structure in Prob. Fields used:
%      IterPrint  Print short information each iteration
%      eps_g     Gradient convergence tolerance
%      eps_x     Convergence tolerance in x
%      size_x    Approximate size of optimal variable values
%      MaxIter   Maximal number of iterations
%      wait      Pause after printout if true
%      cTol      Constraint violation convergence tolerance
%      bTol      Linear constraint violation convergence tolerance
%      xTol      Variable violation tolerance
%      QN_InitMatrix  Initial Quasi-Newton matrix, if not empty,
%                     otherwise use identity matrix
%
% OUTPUT PARAMETERS
% Result      Structure with results from optimization
%    x_k      Optimal point
%    v_k      Lagrange multipliers NOT USED
%    f_k      Function value at optimum
%    g_k      Gradient vector at optimum
%    x_0      Starting value vector
%    f_0      Function value at start
%    c_k      Constraint values at optimum
%    cJac     Constraint derivative values at optimum
%    xState   Variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
%    bState   Linear constraint: Inactive==0; On lower bound == 1;
%             On upper bound == 2; Equality == 3;
%    cState   State of each general constraint.
%    Iter     Number of iterations
%    ExitFlag Flag giving exit status.
%             0 = Convergence. Constraints fulfilled, small step in x.
%             1 = Infeasible problem?
%             2 = Maximal number of iterations reached (Inform == 101)
%             3 = No progress in either function value or constraint reduction
%    ExitTest Text string giving ExitFlag and Inform information
%    Inform   Code telling type of convergence (ExitFlag=0)
%             1   Iteration points are close
%             2   Small search direction
%             3   Function value below given estimate
%                 Restart with lower Prob.f_Low if minimum not reached
%             4   Projected gradient small
%            10   Karush-Kuhn-Tucker conditions fulfilled

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written Jan 26, 1998.   Last modified Sep 11, 2008.

function Result = nlpSolve(Prob, varargin)

if nargin < 1
    error('nlpSolve needs input structure Prob');
end

solvType=checkType('con');

Prob=ProbCheck(Prob,'nlpSolve',solvType);

PriLev    = Prob.PriLevOpt;   % Print level
[Alg, SolvAlg, B_k]=NameAlg(Prob,Prob.FUNCS.H,PriLev);
if Alg == 2 | Alg == 4
    Prob = iniSolve(Prob,solvType,1,2);
else
    Prob = iniSolve(Prob,solvType,2,2);
end

optParam=Prob.optParam;

% Pick up input variables from Prob structure

[n, x_k, x_km1, xEqual, x_L, x_U, Prob] = BoundInit(Prob);

% Linear constraints

[mA, Ax, bEqual, b_L, b_U, A] = LinearConstr(Prob);

% Nonlinear constraints

[m, c_k, dc_k, cEqual, c_L, c_U] = NonlinConstr(Prob, varargin{:});

mTot=mA+m;

eps_g =  optParam.eps_g;    % Gradient convergence tolerance
eps_x =  optParam.eps_x;    % Convergence tolerance in x
MaxIter= optParam.MaxIter;  % Maximal number of iterations
wait =   optParam.wait;     % Pause after printout if true

cTol =   optParam.cTol;    % Constraint violation convergence tolerance
bTol =   optParam.bTol;    % Linear constraint violation convergence tolerance
xTol =   optParam.xTol;    % Variable violation tolerance
fLow =   Prob.f_Low;

zTol = [bTol*ones(mA,1);cTol*ones(m,1)];

IterPrint = optParam.IterPrint;   % Print short information each iteration

if MaxIter <= 0, MaxIter=100*n; end

DEBUG=0;
k = 0;                      % Iteration index

Result=ResultDef(Prob);
Result.Solver='nlpSolve';
Result.SolverAlgorithm=SolvAlg;

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

% Setup Structure used in QP calls
ProbQP = CreateProbQP(Prob, 'qp', max(500,3*n), PriLev-3, optParam.wait);

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

% Setup Structure used in LP calls
ProbFP = CreateProbQP(Prob, 'fp', max(1000,3*n), PriLev-3, optParam.wait);

Prob.nState = 1;
Prob.Mode   = 2;
f_k = nlp_f(x_k, Prob, varargin{:});      % Function value
Prob.Mode   = 1;
g_k = nlp_g(x_k, Prob, varargin{:});      % Gradient

Result.x_0=x_k;
Result.f_0=f_k;

if Alg==1 | Alg==3
    W_k = nlp_H(x_k, Prob, varargin{:});   % Hessian
    if isempty(W_k)
        W_k=eye(n);
        Alg=2;
        Result.SolverAlgorithm=[SolvAlg '!!! BFGS used'];
    else
        if all(W_k(:) < eps), W_k=eye(n); end
    end
else
    W_k = B_k;
end
Prob.nState = 0;

bl=cat(1,b_L,c_L);
bu=cat(1,b_U,c_U);

z=cat(1,Ax,c_k);

if length(bl) ~= length(z)
    fprintf('nlpSolve: Misspecified constraints\n');
    fprintf('Length from bounds is %d. ',length(bl));
    fprintf('Actual length %d.\n',length(z));
    error('nlpSolve: Error in input');
end

zEqual=cat(1,bEqual,cEqual);
%iE=cat(1,xEqual,bEqual,cEqual);
%I=find(~iE);
%E=find( iE);

[h_k cErr] = ConstraintError(z-bu,bl-z,Alg);

epsSOC = 0.000001;                   % Tolerance. We can use tols in optParam
rho = 5;                             % Trust-region radius
M = 1000;                            % h_k upper bound parameter
u = M*max(1,h_k);                    % h_k upper bound (Section 3.2, page 8)
FILTER = [-inf u 1 1000000];         % [f_k h_k delta_q my] (Sec. 3.5, p. 13)
RestoreFlag = 0;                     % Restoration phase flag
delta_q = 0;                         % HOW SHOULD THIS BE INIT (page 13)

my = 1E-6;

f_kOld=NaN;
h_kOld=NaN;
rhoOld=NaN;
KT=Inf;

global n_f
if IterPrint
    fprintf('Iter Func            f(x)             |step|     ');
    fprintf('TrustRegion    |gradient|   Constraint   Kuhn-Tucker\n');
    fprintf('     cnt                                           ');
    fprintf('   step         (Projected)   Violations      Error\n');
end
d_k =zeros(n,1);

% =========================================================================
while 1 % Algorithm 3: Filter SQP (Section 3.7, page 15)
    % Display
    if PriLev > 1
        fprintf('==================================================\n');
        fprintf('Iteration no: %4.0f  Function value %30.20f\n',k,f_k);
        if mTot > 0
            fprintf('Sum of infeasibilities             %30.20f\n',...
                norm(max(0,bl-z),1) + norm(max(0,z-bu),1));
        end
        if ~(isinf(KT) | isnan(KT))
            fprintf('Kuhn-Tucker error                  %30.20f\n',KT);
        end
    end
    if PriLev > 2
        xprint(x_k,'x_k:');
        xprinte(g_k,'g_k:');
        xprinte(c_k,'c_k:');
    end
    if wait & PriLev > 1, pause; end

    % Solve problem QP (page 4). This QP might be infeasible.
    %
    %     min   0.5 * x' * F * x + c' * x.  x in R^n
    %       x
    %      s/t   x_L <=   x  <= x_U
    %            b_L <= A x  <= b_U

    ProbQP.QP.F=W_k;
    ProbQP.QP.c=g_k;
    ProbQP.x_L=max(-min(rho)*ones(n,1),x_L-x_k);
    ProbQP.x_U=min(min(rho)*ones(n,1),x_U-x_k);
    ProbQP.x_0=zeros(n,1);
    if mA > 0, Ax=A*x_k; end
    z=cat(1,Ax,c_k);

    if isempty(A)
        ProbQP.A=dc_k;
    elseif isempty(dc_k)
        ProbQP.A=A;
    else
        ProbQP.A=[A;dc_k];
    end
    ProbQP.mLin = size(ProbQP.A,1);
    ProbQP.b_L=bl-z;
    ProbQP.b_U=bu-z;

    if DEBUG
        disp(ProbQP.QP.F)
        disp(ProbQP.QP.c)
        disp([ProbQP.x_L ProbQP.x_0 ProbQP.x_U])
        disp([ProbQP.b_L ProbQP.A ProbQP.b_U])
        pause
    end

    if DEBUG
        disp('Solve basic QP')
    end
    ResultQP = tomRunMini(SolverQP,ProbQP);
    ExitFlag=ResultQP.ExitFlag; % ==2 if Phase I simplex failed

    if DEBUG
        d_kp1SOL=ResultQP.x_k;
        l_kp1SOL=ResultQP.v_k;
        xprint([ExitFlag;d_kp1SOL],'ExitFlag,d_kp1SOL')
        xprint(l_kp1SOL,'l_kp1SOL')
        %pause
    end

    if ExitFlag > 0      % Problem QP is infeasible
        % Find new point x_k in restoration phase (Section 3.3, page 8)
        x_restor_start = x_k;
        if DEBUG
            disp('CALLING qpPhase1')
            xprint(x_k,'x:')
        end
        [ResultP1, d_0, rho] = qpPhase1(x_k,Prob,ProbQP,ProbFP,rho,varargin{:});
        x_km1 = x_k;
        g_km1 = g_k;
        x_k = ResultP1.x_k;
        if DEBUG
            disp('after CALLING qpPhase1')
            xprint(x_k,'x:')
        end

        if ResultP1.ExitFlag ~=0
            if DEBUG
                disp('nlpSolve: Restoration phase in qpPhase1 FAILED. Error flag:')
                disp(ResultP1.ExitFlag)
                disp(d_0)
                disp(rho)
            end
        end
        if DEBUG
            disp('After QP (and possibly qphase1 if infeasible')
            xprint(x_k,'x:')
        end
        %if ResultP1.ExitFlag ==2
        if ResultP1.ExitFlag > 0
            %'disturb point'
            r1 = rand(n,1);
            ix = find(r1 >= 0.5 & isinf(x_U));
            d_k(ix) = 0.001*(r1(ix)-0.5);
            ix = find(r1 >= 0.5 & ~isinf(x_U));
            d_k(ix) = 0.001*(r1(ix)-0.5).*(x_U(ix)-x_restor_start(ix));
            %x_k=max(x_L,min(x_U,x_restor_start+rand(n,1)));
            ix = find(r1 < 0.5 & isinf(x_L));
            d_k(ix) = 0.001*(r1(ix)-0.5);
            ix = find(r1 < 0.5 & ~isinf(x_L));
            d_k(ix) = 0.001*(r1(ix)-0.5).*(x_restor_start(ix)-x_L(ix));
            x_k = x_restor_start + d_k;
            % Adjust solution inside lower and upper variable bound
            x_k = max(x_L,min(x_U,x_k));
            if DEBUG
                disp('Exact zero solution in qpPhase1. Could be rank problems.');
                disp('Disturb starting point with random values in [0 1]');
                disp(x_restor_start)
                disp(x_k)
            end
        else
            d_k = x_k - x_restor_start;
        end
        %d_k = x_k - x_restor_start;

        Prob.Mode = 2;
        f_k = nlp_f(x_k, Prob, varargin{:});
        c_k = nlp_c(x_k, Prob, varargin{:});
        if mA > 0, Ax=A*x_k; end
        z=cat(1,Ax,c_k);
        [h_k cErr] = ConstraintError(z-bu,bl-z,Alg);
        Prob.Mode   = 1;
        g_k = nlp_g(x_k, Prob, varargin{:});
        dc_k = nlp_dc(x_k, Prob, varargin{:});
        if Alg==1 | Alg==3
            W_k = nlp_H(x_k,Prob, varargin{:});
            if all(W_k(:) < eps), W_k=eye(n); end
        else
            W_k = updateW(W_k,x_k-x_km1,g_k-g_km1,k,PriLev,eps_x,eps_g);
        end
        RestoreFlag = 1;  % Yes, that's true, a restoration phase has been done

    else  % QP (page 4) solved
        % Solution from QP
        d_kp1=ResultQP.x_k;
        l_kp1=ResultQP.v_k;

        % Compute data for k+1: - Function value f_kp1
        %                       - L1-norm of constr violation h_kp1 (page 3).
        x_kp1 = x_k + d_kp1;
        % Adjust solution inside lower and upper variable bound
        x_kp1 = max(x_L,min(x_U,x_kp1));
        Prob.Mode = 2;
        f_kp1 = nlp_f(x_kp1, Prob, varargin{:});
        c_kp1 = nlp_c(x_kp1, Prob, varargin{:});
        if mA > 0, Ax=A*x_kp1; end
        z=cat(1,Ax,c_kp1);
        [h_kp1 cErrkp1] = ConstraintError(z-bu,bl-z,Alg);

        delta_q=-(g_k'*d_kp1+ d_kp1'*W_k*d_kp1); % (page 13)

        % Check if (f_kp1,h_kp1) is acceptable to the filter
        if DEBUG
            %FILTER
            disp(f_kp1)
            disp(h_kp1)
        end
        if Alg < 3
            accept = CheckFilter(FILTER, f_kp1, h_kp1, cTol);
        else
            accept = CheckFilter1(FILTER, f_kp1, h_kp1, cTol);
        end
        %pause

        if accept % (f_kp1,h_kp1) is acceptable to the filter
            % Accept x_kp1 as new point
            if DEBUG
                disp('Accept x_kp1 as new point')
                xprint([f_k;x_k;h_k],'old f,x,h')
                xprint([f_kp1;x_kp1;h_kp1],'New f,x,h')
                pause
            end
            x_km1 = x_k;
            g_km1 = g_k;
            %
            d_k = d_kp1;
            x_k = x_kp1;
            l_k = l_kp1;                 % IS THIS OK??? (page 5)
            ActRed= f_k - f_kp1;
            delta_q=-(g_k'*d_k+ d_k'*W_k*d_k); % (page 13)
            if mTot > 0
                l_kinf=norm(l_k(n+1:n+mTot),inf);
            else
                l_kinf=0;
            end
            if l_kinf < 1E-6
                my = 1E-6;
            else
                my = max(1E-6,min(1E6,10^(ceil(log10(l_kinf))))); % (page 13)
            end
            f_k = f_kp1;
            c_k = c_kp1;
            h_k = h_kp1;
            Prob.Mode = 1;
            g_k = nlp_g(x_k, Prob, varargin{:});
            dc_k = nlp_dc(x_k, Prob, varargin{:});
            if Alg==1 | Alg==3
                W_k = nlp_H(x_k, Prob, varargin{:});
                if all(W_k(:) < eps), W_k=eye(n); end
            else
                W_k = updateW(W_k,x_k-x_km1,g_k-g_km1,k,PriLev,eps_x,eps_g);
            end

            % Check if we should Decrease or not
            %if delta_q >= delta*h_k^2  + extra
            %end

            % Add (f_kp1,h_kp1) to filter and remove any
            % points that are dominated by (f_kp1,h_kp1)
            if Alg < 3
                FILTER = ChangeFilter(FILTER, f_k, h_k, delta_q, my, ActRed);
            else
                FILTER = ChangeFilter1(FILTER, f_k, h_k, delta_q, my);
            end

            % Possibly increase the trust-region radius (item 2, page 16)
            if abs(rho-norm(d_k,inf)) < xTol
                rho = 2*rho;
            end

        else  % (f_kp1,h_kp1) is NOT acceptable to the filter
            if DEBUG
                disp('NOT ACCEPT POINT')
                disp(f_kp1)
                disp(h_kp1)
                disp('DO 2nd order correction')
            end
            if DEBUG
                disp('CALLING SOC')
                xprint(x_k,'x:')
            end
            % Solve a sequence of QPs for a Second Order Correction (SOC)
            % step, d_hat, and set x_hat=x_k+d_hat (Section 3.1, page 7)
            %[x_hat,f_hat,g_hat,W_hat,c_hat,dc_hat,h_hat,l_hat,r, ...
            %  accept,SOCerr] = SOC(x_k,g_k,Alg,W_k,dc_k,d_kp1, ...
            %                   bl,bu,c_L,c_U,b_L,b_U,A,x_L,x_U,FILTER,rho, ...
            %                   epsSOC,ProbQP,Prob,varargin{:});
            [x_hat,r,accept,SOCerr] = SOC(x_k,g_k,Alg,W_k,dc_k,d_kp1, ...
                bl,bu,c_L,c_U,b_L,b_U,A,x_L,x_U,FILTER,rho, ...
                epsSOC,SolverQP,ProbQP,Prob,varargin{:});

            if accept % (f_hat,h_hat) is acceptable to the filter
                if DEBUG
                    disp('SOC point is accepted by the filter')
                end
                % Accept x_hat as new point
                x_km1 = x_k;
                if DEBUG
                    disp('ACCEPT afterCALLING SOC')
                    xprint(x_k,'x:')
                end
                g_km1 = g_k;
                %
                d_k = x_hat - x_k;
                if DEBUG
                    xprint(d_k,'d_k')
                end
                x_k = x_hat;
                % Do not update l_k (OK?)
                % Do not update delta_q and my (OK?)

                % SEPT: HKH NOW ADDED THESE TWO LINES
                Prob.Mode = 2;
                f_k = nlp_f(x_k, Prob, varargin{:});
                c_k = nlp_c(x_k, Prob, varargin{:});

                if mA > 0, Ax=A*x_k; end
                z=cat(1,Ax,c_k);
                [h_k cErr] = ConstraintError(z-bu,bl-z,Alg);
                Prob.Mode = 1;
                g_k = nlp_g(x_k, Prob, varargin{:});
                dc_k = nlp_dc(x_k, Prob, varargin{:});
                if Alg==1 | Alg==3
                    W_k = nlp_H(x_k, Prob, varargin{:});
                    if all(W_k(:) < eps), W_k=eye(n); end
                else
                    W_k = updateW(W_k,x_k-x_km1,g_k-g_km1,k,PriLev,eps_x,eps_g);
                end

                %f_k = f_hat;
                %c_k = c_hat;
                %h_k = h_hat;
                %dc_k = dc_hat;
                %g_k = g_hat;
                %if Alg==1 | Alg==3
                %   W_k = W_hat;  % Computed in SOC
                %else
                %   W_k = updateW(W_k,x_k-x_km1,g_k-g_km1,k,PriLev,eps_x,eps_g);
                %end

                % Add (f_hat,h_hat) to filter and remove any
                % points that are dominated by (f_hat,h_hat)
                if Alg < 3
                    FILTER = ChangeFilter(FILTER, f_k, h_k, delta_q, my, ActRed);
                else
                    FILTER = ChangeFilter1(FILTER, f_k, h_k, delta_q, my);
                end

                % Possibly increase trust-region radius (item 2, page 16)
                if r <= 0.1
                    rho = 2*rho;
                    if DEBUG
                        disp('Increase rho')
                        disp(rho)
                    end
                end

            else
                if RestoreFlag % First iteration after restoration phase
                    if DEBUG
                        disp('% Accept the best SOC step (set x_k=x_hat)')
                        xprinte(x_k,'x_k:')
                    end

                    % Accept the best SOC step (set x_k=x_hat)
                    x_km1 = x_k;
                    if DEBUG
                        disp('1st iter afterCALLING SOC')
                        xprint(x_k,'x:')
                    end
                    g_km1 = g_k;
                    %
                    d_k = x_hat - x_k;
                    x_k = x_hat;
                    % Do not update l_k (OK?)
                    % Do not update delta_q and my (OK?)
                    % SEPT: HKH NOW ADDED THESE LINES
                    Prob.Mode = 2;
                    f_k = nlp_f(x_k, Prob, varargin{:});
                    c_k = nlp_c(x_k, Prob, varargin{:});
                    if mA > 0, Ax=A*x_k; end
                    z=cat(1,Ax,c_k);
                    [h_k cErr] = ConstraintError(z-bu,bl-z,Alg);
                    Prob.Mode = 1;
                    g_k = nlp_g(x_k, Prob, varargin{:});
                    dc_k = nlp_dc(x_k, Prob, varargin{:});
                    if Alg==1 | Alg==3
                        W_k = nlp_H(x_k, Prob, varargin{:});
                        if all(W_k(:) < eps), W_k=eye(n); end
                    else
                        W_k = updateW(W_k,x_k-x_km1,g_k-g_km1,k,PriLev,eps_x,eps_g);
                    end

                    %f_k = f_hat;
                    %g_k = g_hat;
                    %c_k = c_hat;
                    %h_k = h_hat;
                    %dc_k = dc_hat;
                    %if Alg==1 | Alg==3
                    %   W_k = W_hat;  % Computed in SOC
                    %else
                    %   W_k = updateW(W_k,x_k-x_km1,g_k-g_km1,k,PriLev,eps_x,eps_g);
                    %end

                    if h_k > u
                        disp('BAD RESTORATION')
                        disp(cErr)
                        disp(h_k)
                        disp(u)
                    else
                        % Remove all blocking entries that dominates
                        % (f_hat,g_hat) from the filter (Section 3.4)
                        %FILTER = FILTER(find(FILTER(:,1) + xTol > f_k | ...
                        %         FILTER(:,2) + cTol > h_k),:);
                        FILTER = FILTER(find(FILTER(:,1)  > f_k | ...
                            FILTER(:,2) > h_k),:);
                        % Reduce the upper bound u, add (-inf,u) to filter
                        % and and remove points dominated by (-inf,u)

                        u = max(h_k,u/10);
                    end

                    % Add (f_hat,h_hat) to filter? No, not if bellow works


                    if Alg < 3
                        FILTER = ChangeFilter(FILTER, -inf, u, 1, 1000000, ActRed);
                    else
                        FILTER = ChangeFilter1(FILTER, -inf, u, 1, 1000000);
                    end

                else
                    if DEBUG
                        disp('% Reject the  SOC step')
                        xprinte(x_k,'x_k:')
                    end
                    % Reject the step (x_kp1=x_k, i.e. do not change x_k)
                    % Reduce the trust region radius (item 1, page 16)
                    conver = 0; % Convergence flag
                    if DEBUG
                        disp('Reject afterCALLING SOC')
                        xprint(x_k,'x:')
                    end

                    % HKH added safety row
                    if isempty(d_k), d_k=d_kp1; end

                    nd_k=norm(d_k,inf);
                    % Safeguard update of rho
                    if nd_k~=0
                        rho = max([1E-1*eps,1E-2*rho,min(rho,nd_k)/2]);
                    else
                        rho = rho/2;
                    end
                    if DEBUG
                        disp(rho)
                    end
                    if rho <= 0.01*eps_x
                        d_k = zeros(n,1);  % Stop meaningless iterations
                        if DEBUG
                            disp('SET d_k to 0')
                        end
                    end

                end
            end
        end
        RestoreFlag = 0;      % Not the restoration phase
    end
    if DEBUG
        disp(rho)
        disp([x_L x_k x_U])
    end

    k = k + 1;
    % Check convergence criterion (item 5, page 16)
    % Maybe not necc. to compute z and cErr here
    z=cat(1,Ax,c_k);
    cErr=max(z-bu,bl-z);

    if DEBUG
        xprinte(g_k,'g_k:')
        xprinte(l_k,'l_k:')
        xprinte(z,'z:  ')
    end

    Result.x_k=x_k;
    Result.g_k=g_k;
    Result.c_k=c_k;
    Result.cJac=dc_k;

    [v_k, Zv, P, Z_L, ccErr, ceq, cineq, gProj] = LagMult(Prob,Result);

    if mTot == 1
        Anorm=norm(ProbQP.A);
    elseif mTot==0
        Anorm=[];
    else
        Anorm=sqrt(sum(ProbQP.A'.^2));
    end
    if mTot > 0
        myMax=max([norm(g_k),max(abs(v_k(1:n))) max(Anorm(:).*v_k(n+1:n+mTot))]);
        if myMax==0, myMax=1; end

        KT=norm(g_k-v_k(1:n)-ProbQP.A'*v_k(n+1:n+mTot))/myMax;
    else
        myMax=norm(g_k);
        if myMax==0, myMax=1; end
        KT=norm((g_k-v_k(1:n))/myMax);
    end
    %fprintf('Kuhn-Tucker error %30.20f\n',max(KT));

    if KT <= cTol & ( mTot==0 | all(cErr <= zTol.*max(1,abs(z))) )
        if DEBUG
            disp([x_L x_k x_U v_k(1:n)])
            disp([bl z bu v_k(n+1:n+mTot)])
        end
        xT=xTol*max(1,abs(x_k));
        iL=abs(x_L-x_k) < xT & ~xEqual;
        iU=abs(x_U-x_k) < xT & ~xEqual;
        iF=find(~(iL | iU | xEqual));
        KKT=all(v_k(iL) >= -xT(iL)) & all(v_k(iU) <= xT(iU)) & ...
            all(abs(v_k(iF)) < xT(iF));
        if KKT & mTot > 0
            zT=zTol.*max(1,abs(z));
            iL=abs(bl-z) < zT & ~zEqual;
            iU=abs(bu-z) < zT & ~zEqual;
            iF=find(~(iL | iU | zEqual));
            KKT=all(v_k(n+find(iL)) > -zT(find(iL))) & ...
                all(v_k(n+find(iU)) <  zT(find(iU))) & ...
                all(abs(v_k(n+iF)) < zT(iF));
        end
    else
        KKT=0;
    end

    %if k > 20
    %pause
    %end
    stop = 0;
    if mA > 0, Ax=A*x_k; end
    %&conver &accept

    if IterPrint
        fprintf('%4d %4d %25.17f %11.7e %10.6f %17.7e %10.5e %10.5e\n', ...
            k, n_f, f_k, norm(d_k), rho, norm(gProj), h_k, KT);
    end
    % HKH
    % fprintf('%30.18f %30.18f\n',norm(d_k),rho);


    if all(cErr <= zTol.*max(1,abs(z)))
        % Item 5 is using eps*sqrt(n+m)
        %if norm(d_k, inf)<=eps_x % MORE !!!

        %if KT <= cTol  &  max(abs(d_k./max(abs(x_k),size_x))) <= eps_x
        if KT <= cTol & KKT
            stop = 1;
            Inform = 10;
            Flag=0;
        elseif max(abs(d_k./max(1,abs(x_k)))) <= eps_x  & ...
                max(abs(gProj./max(1,abs(x_k)))) <= 100*eps_x
            % Item 5 is using eps*sqrt(n+m) and epd*sqrt(n)
            stop = 1;
            Inform = 1;
            Flag=0;
        elseif max(abs(gProj./max(1,abs(x_k)))) <= eps_x & ...
                max(abs(d_k./max(1,abs(x_k)))) <= 100*eps_x
            % Test on the projected gradient
            % Accept if low projected gradient, and almost convergence in x_k
            stop = 1;
            Inform = 4;
            Flag=0;
        elseif f_k <= fLow
            % Below the target function value
            stop = 1;
            Inform = 3;
            Flag=0;
        elseif max(abs(d_k./max(1,abs(x_k)))) <= eps_x
            % HKH
            % Could do some kind of reset here
            % 'Small d_k - No convergence yet'
            %  rho = 5;
            %  d_k = (d_k+1) * 1E15;
        end
    elseif h_k <= cTol
        disp(rho)
        disp('******************* NOT TESTING**********************')
    end
    if ~stop
        if all(abs(d_k) <= xTol) & h_k <= cTol
            stop = 1;
            Inform = 2;
            Flag=0;
        elseif all(abs(d_k) <= xTol) & ~RestoreFlag & 0
            % Infeasible problem?
            stop = 1;
            Inform = 1;
            Flag=1;
        elseif h_k <= cTol & f_k <= fLow
            % Below the target function value
            stop = 1;
            Inform = 3;
            Flag=0;
        elseif abs(h_k-h_kOld) < xTol & abs(f_k-f_kOld) < xTol & ...
                abs(rho-rhoOld) < xTol
            % No progress at all!  STOP
            rho=rho/4;
            %stop = 1;
            %Inform = 3;
            %Flag=3;
        end
    end
    if stop > 0 | k >= MaxIter
        if k>=MaxIter
            Inform=101;
            Flag=2;
        end
        Result.x_k=x_k;
        Result.v_k=v_k;
        Result.Iter=k;
        Result.ExitFlag=Flag;
        Result.Inform=Inform;
        Result.ExitText=ExitText(Flag,Inform);
        Result.f_k=f_k;
        Result.g_k=g_k;
        Result.c_k=c_k;
        Result.cJac=dc_k;
        if Alg==1 | Alg==3
            Result.H_k=W_k;
        else
            Result.B_k=W_k;
        end
        % State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
        Result.xState=double( (x_k==x_L)+2*(x_k==x_U) );
        if ~isempty(A)
            Result.bState=double( (Ax-bTol <= b_L)+2*(Ax >= b_U-bTol) );
        else
            Result.bState=[];
        end
        if ~isempty(c_k)
            Result.cState=double( (c_k <= c_L+cTol)+2*(c_k >= c_U-cTol) );
        else
            Result.cState=[];
        end
        Result=endSolve(Prob,Result);
        return                  % Convergence
    end
    rhoOld=rho;
    %if abs(d_k) < 1E-12 & rho < 1E-8
    %   % HKH
    %   % This is a weird condtion, no convergence, and small trust region
    %   keyboard
    %end
    h_kOld=h_k;
    f_kOld=f_k;
end % Main loop ends ========================================================
% ===========================================================================


% ---------------------------------------------------------------------------
function accept = CheckFilter(FILTER, f_k, h_k, cTol)
% ---------------------------------------------------------------------------

% Check if point (f_k,h_k) is acceptable to the filter ----------------------

accept = 1;                       % Accept flag (default: point accepted)
%plot(FILTER(1:end-1,2),FILTER(1:end-1,1),'*',h_k,f_k,'o')

if isempty(FILTER)
    return                         % Terminate function
end

DEBUG=1;
if DEBUG
    PrintMatrix(FILTER,'FILTER')
    xprint([f_k h_k],'f_k,h_k')
    %pause
end

tol=1E-12;
delta = 1E-6;
beta = 0.99;                      % Eq (3), page 13
alpha1 = 0.01;                    % Eq (4), page 13
alpha2 = 0.00001;                 % Eq (4), page 13

% Check Eq. (3) and eq. (4), page 13 (FILTER = [f_k h_k delta_q my])
if h_k <= cTol
    eq3 = h_k-tol <= 0.1*FILTER(:,2);  % (page 14)
else
    eq3 = h_k-tol <= beta*FILTER(:,2);
end
% New check New (2.4)
eq4=f_k-tol <= FILTER(:,1)-max(alpha1 * FILTER(:,3).* (FILTER(:,3) > ...
    delta*h_k^2), alpha2*FILTER(:,2).*FILTER(:,4));

if ~all(eq3|eq4)
    accept = 0;
    return                         % (3) and (4) failed: terminate function
end

% Check North-West corner rule (Section 3.6, page 14)
if size(FILTER,1)  > 1
    [h_low,index] = min(FILTER(:,2)); % Not needed if filter is ordered
    if h_k+tol <= h_low
        myy = 1000*FILTER(index,4); % Overestimate of penalty parameter
        if f_k+myy*h_k - tol > FILTER(index,1)+myy*FILTER(index,2)
            accept = 0;
            return               % N-W corner rule failed: terminate function
        end
    end
end

% Check South-East corner rule (Section 3.6, page 15)
if size(FILTER,1) > 1
    [f_low,index] = min(isfinite(FILTER(:,1)).*FILTER(:,1));
    if f_k+tol <= f_low
        myy = FILTER(index,4)/1000; % Overestimate of penalty parameter
        if f_k+myy*h_k - tol > FILTER(index,1)+myy*FILTER(index,2)
            accept = 0;
            return               % S-E corner rule failed: terminate function
        end
    end
end

% N-W and S-E violated during unblocking ??? (Section 3.6, page 15)

% Function CheckFilter ends -------------------------------------------------
% ---------------------------------------------------------------------------


% ---------------------------------------------------------------------------
function FILTER = ChangeFilter(FILTER, f_k, h_k, delta_q, my, ActRed)
% ---------------------------------------------------------------------------

%tol=1E-12;
tol=0;

% Accept (f_k,h_k) to the filter and remove points dominated by (f_k,h_k) ---

% Remove points dominated by (f_k,h_k). Definition 1, page 4
FILTER = FILTER(find(f_k - tol >= FILTER(:,1) | h_k - tol >= FILTER(:,2)),:);

% ?????
% WHAT TO SET delta > 0  TO???

delta=1;

fprintf('Delta_g %f h_k %f h_k^2 %f delta %f Skip %d\n', ...
    delta_q,h_k,h_k^2,delta,delta_q >= delta*h_k^2);

% Check if we should update or not
if delta_q >= delta * h_k^2
    return
end

% Accept (f_k,h_k) to the filter
if isempty(FILTER)
    FILTER = [f_k h_k delta_q my];
else
    if f_k + tol >= FILTER(1,1) & h_k - tol <= FILTER(1,2)
        FILTER = [[f_k h_k delta_q my];FILTER];
    elseif f_k - tol <= FILTER(end,1) & h_k + tol >= FILTER(end,2)
        FILTER = [FILTER;[f_k h_k delta_q my]];
    else
        index = find(f_k + tol >= FILTER(:,1) & h_k - tol <= FILTER(:,2));
        FILTER=[FILTER(1:index-1,:);[f_k h_k delta_q my];FILTER(index:end,:)];
    end
end


% Function ChangeFilter ends ------------------------------------------------
% ---------------------------------------------------------------------------


% ---------------------------------------------------------------------------
% Second Order Correction (SOC) step (Section 3.1, page 7) ------------------
%function [x_hat,f_hat,g_hat,W_hat,c_hat,dc_hat,h_hat,l_hat,r,accept,SOCerr] ...
%          = SOC(x_k,g_k,Alg,W_k,dc_k,d_kp1, bl,bu, c_L,c_U,b_L,b_U,A,...
%                x_L,x_U,FILTER,rho,epsSOC,ProbQP,Prob,varargin)

function [x_hat,r,accept,SOCerr] ...
    = SOC(x_k,g_k,Alg,W_k,dc_k,d_kp1, bl,bu, c_L,c_U,b_L,b_U,A,...
    x_L,x_U,FILTER,rho,epsSOC,SolverQP,ProbQP,Prob,varargin)

nargin;
DEBUG=0;

d_hat = d_kp1;
x_hat = x_k + d_hat;
% Adjust solution inside lower and upper variable bound
x_hat = max(x_L,min(x_U,x_hat));
Prob.Mode = 0;
c_hat = nlp_c(x_hat, Prob, varargin{:});
h_hat = inf;

n     = length(x_k);
cTol = Prob.optParam.cTol;
Ax    = zeros(0,1);

% We should run a sequence of SOC steps. But how shall we update? The SOC
% step is described on page 395 in the book "Practical Methods of
% Optimization" by Fletcher.

% These items are NOT changed
ProbQP.QP.F=W_k;
ProbQP.QP.c=g_k;
ProbQP.x_L=max(-min(rho)*ones(n,1),x_L-x_k);
ProbQP.x_U=min( min(rho)*ones(n,1),x_U-x_k);
ProbQP.x_0=zeros(n,1);

Iter = 1;
while 1
    % Solve problem QP2 (page 7)
    % Hot start: d_hat?
    % Solve problem QP
    if isempty(A)
        ProbQP.A=dc_k;
    elseif isempty(dc_k)
        ProbQP.A=A;
    else
        ProbQP.A=[A;dc_k];
    end
    ProbQP.mLin = size(ProbQP.A,1);
    Ad=ProbQP.A*d_hat;

    if size(A,1) > 0, Ax = A*x_hat; end
    z=cat(1,Ax,c_hat(:));

    ProbQP.b_L=bl-z + Ad;
    ProbQP.b_U=bu-z + Ad;


    if DEBUG
        disp('Solve SOC step QP')
    end
    ResultQP = tomRunMini(SolverQP,ProbQP);

    d_hat=ResultQP.x_k;
    SOCerr=ResultQP.ExitFlag;
    DEBUG=0;
    if DEBUG
        disp('nlpSolve, SOC: Result of QP');
        disp(SOCerr)
        disp(Iter)
        xprinte(d_hat,'d_hat:');
        %%pause
    end

    % Criterion 2, page 8
    if SOCerr > 0
        accept = 0;
        r=NaN;
        if Iter == 1
            x_hat = x_k;

            %f_hat = nlp_f(x_hat, Prob, varargin{:});
            %g_hat = g_k;
            %W_hat = W_k;
            %dc_hat = dc_k;
            %c_hat = nlp_c(x_hat, Prob, varargin{:});
            %if size(A,1) > 0, Ax=A*x_hat; end
            %z=cat(1,Ax,c_hat(:));

            %h_hat = norm(max(0,bl-z),1) + norm(max(0,z-bu),1);
            %l_hat = NaN; % NOT USED
        end
        return                   % Terminate function
    end

    % Update parameters (x_hat = accepted point)
    h_hatm1 = h_hat;
    x_hat = x_k + d_hat;
    % Adjust solution inside lower and upper variable bound
    x_hat = max(x_L,min(x_U,x_hat));
    Prob.Mode = 0;
    f_hat = nlp_f(x_hat, Prob, varargin{:});
    c_hat = nlp_c(x_hat, Prob, varargin{:});
    if size(A,1) > 0, Ax=A*x_hat; end
    z=cat(1,Ax,c_hat);

    if ~isempty(z)
        h_hat = norm(max(0,bl-z),1) + norm(max(0,z-bu),1);
    else
        h_hat = 0;
    end

    %g_hat = nlp_g(x_hat, Prob, varargin{:});
    %dc_hat = nlp_dc(x_hat, Prob, varargin{:});
    %if Alg==1 | Alg==3
    %   % Use Hessian in new point
    %   W_hat = nlp_H(x_hat, Prob, varargin{:});
    %else
    %   W_hat = W_k;
    %end

    r = h_hat/h_hatm1;          % SOC convergence rate

    % Check if (f_hat,h_hat) is acceptable to the filter
    if Alg < 3
        accept = CheckFilter(FILTER, f_hat, h_hat, cTol);
    else
        accept = CheckFilter1(FILTER, f_hat, h_hat, cTol);
    end
    % Criterion 1, 3 and 4 on page 8

    if accept| r>=0.25 | h_hat < epsSOC  | Iter > 20
        return                   % Terminate function
    end
    Iter=Iter+1;
end
% Function SOC ends ---------------------------------------------------------
% ---------------------------------------------------------------------------

% =========================================
function [Alg, SolvAlg, B_k]=NameAlg(Prob,userH,PriLev)
% =========================================

% Alg = 1 = Fletcher-Leyffer SQP-Filter using Hessian, numerical or analytical
% Alg = 2 = Fletcher-Leyffer SQP-Filter, BFGS update
% Alg = 3 = Fletcher-Leyffer Filter SQP  using Hessian, numerical or analytical
% Alg = 4 = Fletcher-Leyffer Filter SQP, BFGS update

Alg=Prob.Solver.Alg;
if isempty(Alg), Alg=0; end

% FIX NOW, DO NOT RUN NEW SQP FILTER. OBZ HKH FIX
if Alg == 1, Alg = 3; end
if Alg == 2, Alg = 4; end

if Alg==1 & Prob.NumDiff==0 & (~isempty(userH) | any(Prob.ADObj == -1))
    %SolvAlg='Fletcher-Leyffer SQP-Filter';
    SolvAlg='Fletcher-Leyffer Filter SQP. Analytic Hessian';
elseif Alg==1
    %SolvAlg='Fletcher-Leyffer SQP-Filter. Numerical Hessian';
    SolvAlg='Fletcher-Leyffer Filter SQP. Numerical Hessian';
elseif Alg==2
    %SolvAlg='Fletcher-Leyffer SQP-Filter. BFGS update';
    SolvAlg='Fletcher-Leyffer Filter SQP. BFGS update';
elseif Alg==3 & Prob.NumDiff==0 & (~isempty(userH) | any(Prob.ADObj == -1))
    SolvAlg='Fletcher-Leyffer Filter SQP. Analytic Hessian';
elseif Alg==3
    SolvAlg='Fletcher-Leyffer Filter SQP. Numerical Hessian';
elseif Alg==4
    SolvAlg='Fletcher-Leyffer Filter SQP. BFGS update';
else
    if Prob.NumDiff==0 & (~isempty(userH) | any(Prob.ADObj == -1))
        SolvAlg='Fletcher-Leyffer Filter SQP. Analytic Hessian';
        Alg=3;
    else
        SolvAlg='Fletcher-Leyffer Filter SQP. BFGS update';
        Alg=4;
    end
end

% Find initial Quasi-Newton matrix
n = Prob.N;
if Alg==2 | Alg==4
    if ~isempty(Prob.optParam.QN_InitMatrix)
        B_k = Prob.optParam.QN_InitMatrix;   % B = User given matrix
        if size(B_k,1)~=n | size(B_k,2)~=n
            if PriLev >=0
                disp('Illegal dimensions on initial Quasi Newton matrix');
                disp('ucSolve is using identity matrix instead');
            end

            B_k = eye(n);                % B = Use identity matrix
        end
    else
        B_k = eye(n); % B = Start Quasi-Newton with identity matrix
    end
else
    B_k=[];
end

% ---------------------------------------------------------------------------
% Safeguarded BFGS update of approximate Hessian ----------------------------

function B_k = updateW(B_k,z,y,k,PriLev,eps_x,eps_g)
% z is z=x_k-x_km1; If using idx~=[1:n], then z=z(idx)
% y is y=g_k-g_km1; If using idx~=[1:n], then y=y(idx)

if norm(z) > eps_x & norm(y) > eps_g
    n=length(z);
    idx=1:n;
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
        fprintf('BFGS Hessian update iter %3.0f\n',k);
        fprintf('zw %10.6e zBz %10.6e \n',zw,zBz)
    end
    if zw < 1E-13 | zBz < 1E-13
        if PriLev > -2
            fprintf('BFGS Hessian update step dangerous!\n');
            fprintf('max(w)/zw %20.10e\n',max(w)/zw);
            fprintf('min(w)/zw %20.10e\n',min(w)/zw);
            fprintf('max(Bz)/zBz %20.10e\n',max(Bz)/zBz);
            fprintf('min(Bz)/zBz %20.10e\n',min(Bz)/zBz);
        end
    else
        B_k(idx,idx)=B_k(idx,idx)+(w/zw)*w'-(Bz/zBz)*Bz';
    end
    if PriLev > 3
        fprintf('W_k matrix after BFGS update\n');
        PrintMatrix(B_k,'W_k:')
    end
    if PriLev > 2
        fprintf('W_k matrix eigenvalues:\n');
        xprint(eig(B_k),'eig:')
    end
end

% Function updateW ----------------------------------------------------------
% ---------------------------------------------------------------------------

% ---------------------------------------------------------------------------
% Check if point (f_k,h_k) is acceptable to the filter ----------------------
function accept = CheckFilter1(FILTER, f_k, h_k, cTol)

tol = 1E-12;

accept = 1;                       % Accept flag (default: point accepted)

%plot(FILTER(1:end-1,2),FILTER(1:end-1,1),'*',h_k,f_k,'o')

if isempty(FILTER)
    return                         % Terminate function
end

% disp(FILTER)
% fprintf('%20.15f %20.15f %15.10f %15.10f\n',FILTER);

beta   = 0.99;                    % Eq (3), page 13
alpha1 = 0.01;                    % Eq (4), page 13
alpha2 = 0.00001;                 % Eq (4), page 13

% Check Eq. (3) and eq. (4), page 13 (FILTER = [f_k h_k delta_q my])
if h_k <= cTol
    eq3 = h_k-tol <= 0.1*FILTER(:,2);  % (page 14)
else
    eq3 = h_k-tol <= beta*FILTER(:,2);
end

eq4=f_k-tol <= FILTER(:,1)-max(alpha1*FILTER(:,3),...
    alpha2*FILTER(:,2).*FILTER(:,4));

if ~all(eq3|eq4)
    accept = 0;
    return                         % (3) and (4) failed: terminate function
end

% Check North-West corner rule (Section 3.6, page 14)
if size(FILTER,1)>1
    [h_low,index] = min(FILTER(:,2)); % Not needed if filter is ordered
    if h_k + tol <= h_low
        myy = 1000*FILTER(index,4); % Overestimate of penalty parameter
        if f_k+myy*h_k-tol > FILTER(index,1)+myy*FILTER(index,2)
            accept = 0;
            return               % N-W corner rule failed: terminate function
        end
    end
end

% Check South-East corner rule (Section 3.6, page 15)
if size(FILTER,1)>1
    [f_low,index] = min(isfinite(FILTER(:,1)).*FILTER(:,1));
    if f_k + tol <= f_low
        myy = FILTER(index,4)/1000; % Overestimate of penalty parameter
        if f_k+myy*h_k-tol > FILTER(index,1)+myy*FILTER(index,2)
            accept = 0;
            return               % S-E corner rule failed: terminate function
        end
    end
end

% N-W and S-E violated during unblocking ??? (Section 3.6, page 15)

% Function CheckFilter1 ends -------------------------------------------------
% ---------------------------------------------------------------------------


% ---------------------------------------------------------------------------
% Accept (f_k,h_k) to the filter and remove points dominated by (f_k,h_k) ---
function FILTER = ChangeFilter1(FILTER, f_k, h_k, delta_q, my)

%tol = 1E-12;
%tol1 = 1E-12;
tol = 0;
tol1 = 0;


% fprintf('I  %20.15f %20.15f %15.10f %15.10f\n',FILTER);
% Remove points dominated by (f_k,h_k). Definition 1, page 4

FILTER = FILTER(find(f_k - tol1 >= FILTER(:,1) | h_k - tol1 >= FILTER(:,2)),:);

% Accept (f_k,h_k) to the filter

if isempty(FILTER)
    FILTER = [f_k h_k delta_q my];
else
    if f_k + tol >= FILTER(1,1) & h_k - tol <= FILTER(1,2)
        FILTER = [[f_k h_k delta_q my];FILTER];
    elseif f_k - tol <= FILTER(end,1) & h_k + tol >= FILTER(end,2)
        FILTER = [FILTER;[f_k h_k delta_q my]];
    else
        index = find(f_k + tol >= FILTER(:,1) & h_k - tol <= FILTER(:,2));
        FILTER=[FILTER(1:index-1,:);[f_k h_k delta_q my];FILTER(index:end,:)];
    end
end

% Function ChangeFilter1 ends ------------------------------------------------
% ---------------------------------------------------------------------------

function [h_k,cErr] = ConstraintError(v1,v2,Alg)

if isempty(v1) & isempty(v2)
    cErr=[];
    h_k=0;
    return
end

cErr=max(v1,v2);

if Alg < 3
    h_k = norm(max(0,cErr),Inf);  % Constr. violation, max norm
else
    h_k = norm(max(0,cErr),1);    % Constr. violation, L1 - sum of abs
end

% MODIFICATION LOG:
%
% 980825  hkh  Added part to make nlpSolve work with new TOMLAB, see NEW PART
% 980826  hkh  Added init of d_k and extra check setting d_k=d_kp1 if empty
% 981005  hkh  Fixed input/output for new TOMLAB. Adding call to QPOPT.
% 981005  edr  Algorithm fixed for new input/output format.
% 981013  hkh  Added call to iniSolve and endSolve
% 981015  hkh  Changed convergence test
% 981016  hkh  Delete arg n from updateW. Input x_k-x_km1 and g_k-g_km1.
%              Safeguard BFGS update for LP problems, or near linear problems,
%              i.e. when y small. Also check for small x_k-x_km1 difference.
% 981019  hkh  Changed all(eq3+eq4) to all(eq3|eq4) (same result)
%              Avoid rho being set to 0 if norm(d_k,inf)==0
%              Also safeguard rho, not getting too small too fast
%              Avoid test of m==0 in setting ProbQP
% 981026  hkh  Use optParam.Method to make choice of QPOPT or qpSolve
%              Use field SolverAlgorithm and fix printing levels
% 981028  hkh  Put f_0 = f(x_0) as field in Result
% 981101  hkh  Remove addition of simple bounds as constraints.
%              Corrected and safeguarded formula for my.
%              Changed update of rho. Simplified constraint handling, like
%              using (bl-z,z-bu) in computations of errors
% 981121  hkh  Error in check if second derivative available
% 981124  hkh  Add check no change in rho/f_k/h_k, then reduce rho.
%

%============================================================================
%
% qpPhase1.m:
%
% function [Result, d_0, rho] = qpPhase1(Prob,rho,varargin)
% function [ResultQ1, d_0, rho] = qpPhase1(x_k,ProbQP,ProbFP,rho,varargin)
%
% Utility routine called from nlpSolve. Implements the restoration phase.
%
% PURPOSE: The constraints in a nonlinear programming problem are defined as
%
%           s/t   x_L <=   x  <= x_U
%                 b_L <= A x  <= b_U
%                 c_L <= c(x) <= c_U.
%
%           In each iteration, using a SQP solver, these constraints are
%           linearized at the relevant x_k, giving a QP. If QP is
%           infeasible, this routine is called to compute a new point x_k
%           which gives a (nearly) feasible QP.
%
% rho  Trust region radius
% d_0  Solution
%

function [ResultQ1, d_0, rho] = qpPhase1(x_k,Prob,ProbQP,ProbFP,rho,varargin)

nargin;
DEBUG=0;

if DEBUG
    disp('*** Start qpPhase1.m')
end

SolverQP = ProbQP.Solver.Name;
SolverFP = ProbFP.Solver.Name;

% Pick up input variables from Prob structure
n   = length(x_k);
d_0 = zeros(n,1);
x_L = ProbQP.x_L(:);
x_U = ProbQP.x_U(:);

b_L = ProbQP.b_L(:);
b_U = ProbQP.b_U(:);
A   = ProbQP.A;

c   = Prob.FUNCS.c;
d2c = Prob.FUNCS.d2c;
c_L = Prob.c_L(:);
c_U = Prob.c_U(:);

Alg=Prob.Solver.Alg;

optParam=Prob.optParam;

ResultQ1=ResultDef(Prob);
ResultQ1.Solver='qpPhase1';

mA=size(A,1);

if mA > 0
    Ax=A*x_k;
    bEqual=b_L==b_U;
else
    Ax=zeros(0,1); % Because b_L(:) makes empty vector be size (0,1)
    bEqual=[];
end
m=max(length(c_L),length(c_U));
N=n+mA+m;

if m > 0
    Prob.Mode = 2;
    c_k  = nlp_c(x_k, Prob, varargin{:});  % Constraints
    Prob.Mode = 1;
    dc_k = nlp_dc(x_k, Prob, varargin{:});% Constraint Jacobian
else
    dc_k = zeros(0,1);
    c_k  = zeros(0,1); % Because c_L(:),c_U(:) makes empty vector be size (0,1)
end

cEqual=c_L==c_U;

bl=[b_L;c_L];
bu=[b_U;c_U];
bcEqual=[bEqual;cEqual];

z=[Ax;c_k];

[h_k cErr] = ConstraintError(z-bu,bl-z,Alg);

PriLev= Prob.PriLevOpt;       % Print level
xTol  =  optParam.xTol;
bTol  =  optParam.bTol;
cTol  =  optParam.cTol;

zTol=[bTol*ones(mA,1);cTol*ones(m,1)];

M       = 1000;                      % h_k upper bound parameter
u       = M*max(1,h_k);              % h_k upper bound (Section 3.2, page 8)
FILTER  = [-inf u 1 1000000];        % [f_k h_k delta_q my] (Sec. 3.5,p. 13)
delta_q = 0;                         % HOW SHOULD THIS BE INIT (page 13)

Jold=zeros(length(bl),1); d_k = []; Jsum=0;

% Define Prob structure for FP - Phase 1 LP problem

ProbFP.x_0   = [];
% Must define a 0-vector of length n c, to use in lp_f to compute c'*x
% ProbFP.QP.c  = [];
ProbFP.QP.c  = zeros(n,1);

% Define Prob structure for Restoration QP problems

%ProbQP      = Prob;
%ProbQP.c_L= [];
%ProbQP.c_U= [];
%ProbQP.p_c= [];
%ProbQP.PriLevOpt=Prob.PriLevOpt-2;

Iter = 0;                               % Iteration index

while rho >= xTol * max(max(1,abs(x_k))) & Iter < 20
    if DEBUG
        fprintf('qpPhase1: Iteration %d. rho %15.9e\n',Iter,rho);
    end

    SOCerr = 0; % SOC error flag

    % Check if problem (page 9) is feasible using FP -  Phase I Simplex.
    %
    %     min   0
    %      p
    %      s/t   max(-rho, x_L - x_k)  <=   p       <= min(rho,x_U - x_k)
    %            b_L -A*x_k            <=   A p     <= b_U - A*x_k
    %            c_L -c_k              <= dc(x_k) p <= c_U - c_k

    % Bounds around x_k are Prob.x_L, Prob.x_U
    ProbFP.x_L  = max(-min(rho)*ones(n,1),Prob.x_L-x_k);
    ProbFP.x_U  = min( min(rho)*ones(n,1),Prob.x_U-x_k);

    if isempty(A)
        if isempty(dc_k)
            ProbFP.A = zeros(0,1);
        else
            ProbFP.A = dc_k;
        end
    elseif isempty(dc_k)
        ProbFP.A    = A;
    else
        ProbFP.A    = [A;dc_k];
    end
    ProbFP.mLin = size(ProbFP.A,1);
    if Iter > 0
        if mA > 0, Ax=A*x_k; end
        z=[Ax;c_k];
    end
    ProbFP.b_L  = bl-z;
    ProbFP.b_U  = bu-z;

    if DEBUG
        disp('qpPhase1: ------------ Run FP - Phase 1 LP --------------')
        disp(rho)
    end

    if rho < 1E-8 & DEBUG
        ProbFP.PriLevOpt=5;
        ProbFP.optParam.wait=1;
    end
    ResultFP = tomRunMini(SolverFP,ProbFP);

    ExitFlag = ResultFP.ExitFlag;
    d_0      = ResultFP.x_k;
    r_k      = ResultFP.r_k;
    v_k      = ResultFP.v_k;
    y        = ResultFP.y_k;
    if DEBUG
        fprintf('LP Phase1 solution. ExitFlag = %d\n',ExitFlag);
        xprint(d_0,'d_0')
        xprint(x_k,'x_k')
        xprint(r_k,'r_k')
        xprint(z,'z  ')
        disp('[x_L d_0 x_U resid v_k]')
        if isempty(r_k), r_k=zeros(N,1); end
        [ProbFP.x_L d_0 ProbFP.x_U r_k(1:n) v_k(1:n)]
        disp('[b_L A*d_0 b_U resid v_k bl z bu]')
        [ProbFP.b_L ProbFP.A*d_0 ProbFP.b_U r_k(n+1:N) v_k(n+1:N) bl z bu]
        pause
    end


    if ExitFlag==0  % LP is feasible
        ResultQ1.x_k = x_k;
        ResultQ1.ExitFlag = 0;
        ResultQ1.Prob = Prob;
        if DEBUG
            disp('*** Leave qpPhase1.m')%DEBUG
        end
        return

    else  % Problem LP (page 9) is not feasible or d_0 nonzero

        %if mA > 0
        %   Ax=A*d_0;
        %   z(1:mA)=Ax;
        %   cErr=max(z-bu,bl-z);
        %end
        %format long
        %disp('resid - computed error')
        %sum(abs(cErr-r_k(n+1:N)))

        % Define set J
        %J=zeros(length(bl),1);
        if isempty(r_k)
            r_k=zeros(N,1);
        end

        % Upper bounds violated are now set to 1
        J= r_k(n+1:N) > zTol .* max(1,abs(z));
        % Lower bounds violated should have -1
        ix = ProbFP.b_L-ProbFP.A*d_0 > zTol .* max(1,abs(z));
        J(ix)=-J(ix);

        % Check on constraint violation
        %J(bu(1+n) + xTol < z(1+n))= 1;   % Upper bounds does not hold
        %J(bl(1+n) - xTol > z(1+n))=-1;   % Lower bounds does not hold
        iJ1=find(J==1);
        iJ2=find(J==-1);
        if isempty(iJ1)
            if isempty(iJ2)
                g_QP=zeros(n,1);
            else
                g_QP=-sum(ProbFP.A(iJ2,:),1);
            end
        elseif isempty(iJ2)
            g_QP=sum(ProbFP.A(iJ1,:),1);
        else
            g_QP=sum(ProbFP.A(find(iJ1),:),1)-sum(ProbFP.A(find(iJ2),:),1);
        end

        if DEBUG
            mPrint(ProbFP.A,'A;dc:')
            xprinti(J,'J:')
            xprinte(g_QP,'g_QP: ')
        end

        % Set the constant vector and solve LP1 to get Lagrange multipliers
        ProbQP.QP.F = [];
        ProbQP.QP.c = full(g_QP(:));
        ProbQP.A    = ProbFP.A;
        ProbQP.b_L  = ProbFP.b_L;
        ProbQP.b_U  = ProbFP.b_U;
        ProbQP.mLin = size(ProbQP.A,1);
        ProbQP.x_L  = ProbFP.x_L;
        ProbQP.x_U  = ProbFP.x_U;
        if rho < 1E-8 & DEBUG
            ProbQP.PriLevOpt=5;
            ProbQP.optParam.wait=1;
        end
        % Change constraint boundaries to Inf where constraint not feasible.
        ProbQP.b_L(iJ2) = -Inf;    % Lower bounds does not hold
        ProbQP.b_U(iJ1) =  Inf;    % Upper bounds does not hold
        % Try a half warm start (could also send the QR and B)
        ProbQP.x_0  = ResultFP.x_k;
        
        if DEBUG
            disp('Solve LP1 in qpPhase')
        end
        ResultQP = tomRunMini(SolverQP,ProbQP);

        B = ResultQP.QP.B;
        d_LP1=ResultQP.x_k;
        l_LP1=ResultQP.v_k;
        y=ResultQP.y_k;

        if DEBUG
            disp('Solution LP1')
            xprinte(l_LP1,'l_LP1: ')
            xprinte(d_LP1,'d_LP1: ')
            pause
        end
        if ResultQP.ExitFlag ~=0 & PriLev > 0
            if DEBUG
                disp('Trouble IN qpPhase, when solving LP1');
                fprintf('ExitFlag = %d\n',ResultQP.ExitFlag);
                rho
            end
            %ExitFlag=ResultQP.ExitFlag;
        end

        % Expanded LP1 if second derivatives exists
        %if m > 0
        %   anyCons=any(r_k(n+mA+1:N)~=0);
        %else
        %   anyCons=0;
        %end
        %if ~anyCons

        if m==0
            d_kp1=d_LP1;
            l_kp1=l_LP1;
        else
            % Build QP (expanded LP1)
            %if 0
            %   W_QP = zeros(n,n);
            %   cons = length(resid_exp);
            %   lambda_QP = ones(cons,1);
            %   for i=1:length(T_exp)
            %      lambda_QP(T_exp(i)) = l_LP1(i);
            %   end
            %   %
            %   for i=1:n
            %      lam_temp = zeros(n,1);
            %      lam_temp(i) = lambda_QP(i);
            %      W_QP = W_QP + nlp_d2c(x_k, lam_temp, Prob, varargin{:});
            %   end
            %   for i=1:n
            %      lam_temp = zeros(n,1);
            %      lam_temp(i) = lambda_QP(n+i);
            %      W_QP = W_QP - nlp_d2c(x_k, lam_temp, Prob, varargin{:});
            %   end
            %end
            % Now use only nonlinear constraints, linear have 2nd der zero

            l_QP=l_LP1(n+mA+1:n+mA+m); % OK for constraints on upper bound
            % Find constraints active on lower bounds
            ix=find(abs(c_k-c_L) < cTol & c_L~=c_U);
            if DEBUG
                if ~isempty(ix)
                    disp('Change signs for Lagrange mults for lower bound constraints')
                    ix
                    pause
                end
            end
            l_QP(ix)=-l_QP(ix); % Change signs on lower bound Lagrange multipliers
            Jc=J(mA+1:mA+m);
            l_QP(Jc== 1)= 1;     % Upper bounds does not hold
            l_QP(Jc==-1)=-1;     % Lower bounds does not hold
            ProbQP.QP.F=nlp_d2c(x_k, l_QP, Prob, varargin{:});
            d2c=ProbQP.QP.F;
            if DEBUG
                xprinti(Jc,'Jc:');
                l_QP
                W_QP=ProbQP.QP.F
                %disp('bl z bu')
                %[bl z bu]
                %PrintMatrix(ProbQP.A,'A;dc:')
                rankW_QP=rank(W_QP)
                lambda=eig(W_QP)
                pause
            end
            %eigV=eig(ProbQP.QP.F);
            %if any(eigV < 0)
            %   disp('qpPhase1: Negative eigenvalues in W_QP');
            %   eigV
            %end

            % Solve problem QP (LP1 with W_kJT on page 9)
            ProbQP.x_0=d_LP1;
            if DEBUG
                disp('Solve Expanded QP-LP1 in qpPhase')
            end
            ResultQP = tomRunMini(SolverQP,ProbQP);

            d_kp1 = ResultQP.x_k;
            l_kp1 = ResultQP.v_k;

            if DEBUG
                disp('Expanded QP solution')
                ExitFlag=ResultQP.ExitFlag;
                xprinte(l_kp1,'l_kp1: ')
                xprinte(d_kp1,'d_kp1: ')
                if ExitFlag > 0
                    pause
                end
            end
        end
        if ResultQP.ExitFlag ~=0 & PriLev > 0

            %disp('qpPhase1: 2nd der QP problem not possible to solve ?????');
            %fprintf('ExitFlag = %d\n',ResultQP.ExitFlag);
            if DEBUG
                %ExitFlag=ResultQP.ExitFlag;
                disp('qpPhase1: 2nd der QP problem not possible to solve ?????');
                fprintf('ExitFlag = %d\n',ResultQP.ExitFlag);
                rho
                pause
            end
        end
        if norm(d_kp1,inf) <= xTol
            d_0=d_kp1;
            ResultQ1.x_k = x_k;
            ResultQ1.ExitFlag = 2;
            if PriLev > 1
                fprintf('qpPhase1: Zero search direction. Failure!');
                fprintf('rho %15.7e Iter %d\n',rho, Iter);
            end
            return
        end

        % Compute data for Iter+1
        x_kp1 = x_k + d_kp1;
        % Adjust solution inside lower and upper variable bound
        x_kp1 = max(Prob.x_L,min(Prob.x_U,x_kp1));
        if mA > 0, Ax=A*x_kp1; end
        % Expand c_k and dc_k giving A_exp'*d <= c_exp
        if m == 0
            c_kp1=zeros(0,1);
        else
            Prob.Mode = 0;
            c_kp1 = nlp_c(x_kp1, Prob, varargin{:});
        end
        z=[Ax;c_kp1];
        cErr=max(z-bu,bl-z);

        if DEBUG
            disp('bl z bu cErr in New Point')
            disp([bl z bu cErr])
            pause
        end

        J0 = J==0;

        if sum(~J0) > 0
            hJ_k   = norm(max(0,cErr(~J0)),1); % Just to update ActRed
            hJ_kp1 = hJ_k;
        else
            hJ_k   = 0;
            hJ_kp1 = 0;
        end
        if sum(J0) > 0
            hT_kp1 = norm(max(0,cErr(J0)),1);
        else
            hT_kp1 = 0;
        end

        % Check if (hJ_kp1,hT_kp1) is acceptable to the filter
        accept = CheckFilter1(FILTER, hJ_kp1, hT_kp1, cTol);
        if DEBUG
            accept
            xprinte(d_kp1,'d_kp1: ')
            if ~accept
                pause
            end
        end

        if accept & ~all(abs(d_kp1) < xTol)
            % (hJ_kp1,hT_kp1) is acceptable to the filter
            % Accept x_kp1 as new point
            d_k = d_kp1;
            x_k = x_kp1;
            l_k = l_kp1;
            if ~isempty(ProbQP.QP.F)
                delta_q=-(g_QP(:)'*d_k+ d_k'*ProbQP.QP.F*d_k); % (page 13)
            else
                delta_q=-g_QP(:)'*d_k; % (page 13)
            end

            l_kinf=norm(l_k,inf);
            if l_kinf < 1E-6
                my = 1E-6;
            else
                my = max(1E-6,min(1E6,10^(ceil(log10(l_kinf))))); % (page 13)
            end

            c_k = c_kp1;
            if m==0
                dc_k=zeros(0,1);
            else
                Prob.Mode = 1;
                dc_k = nlp_dc(x_k, Prob, varargin{:});
            end
            hJ_k = hJ_kp1;
            hT_k = hT_kp1;

            % Add (hJ_kp1,hT_kp1) to filter and remove any
            % points that are dominated by (hJ_kp1,hT_kp1)
            FILTER = ChangeFilter1(FILTER, hJ_k, hT_k, delta_q, my);

            % Possibly increase the trust-region radius (item 2, page 16)
            if abs(rho-norm(d_k,inf)) < xTol
                rho = 2*rho;
            end

        else  % (hJ_kp1,hT_kp1) is NOT acceptable to the filter
            % Solve a sequence of QPs for SOC step
            % [x_hat,c_hat,dc_hat,hJ_hat,hT_hat,r,accept,SOCerr] = ...
            %       SOC2(x_k,dc_k,d_kp1,J_exp,T_exp,W_QP,g_QP,A_QP, ...
            %         c,dc,d2c,bl,bu,c_L,c_U,b_L,b_U,x_L,x_U,FILTER, ...
            %                   rho,eps,ProbQP,Prob,varargin{:});
            %function [x_hat,c_hat,dc_hat,hJ_hat,hT_hat,r,accept,SOCerr] = ...
            %                SOC2(x_k,dc_k,d_kp1,J_exp,T_exp,W_QP,g_QP,A_QP, ...
            %                   c,dc,d2c,bl,bu,c_L,c_U,b_L,b_U,x_L,x_U,FILTER, ...
            %                              rho,eps,ProbQP,Prob,varargin)

            % SOC 2 part -------------------------------------

            d_hat = d_kp1;
            x_hat = x_k + d_hat;
            % Adjust solution inside lower and upper variable bound
            x_hat = max(Prob.x_L,min(Prob.x_U,x_hat));
            if mA > 0, Ax=A*x_hat; end
            Prob.Mode = 0;
            if m==0
                c_hat=zeros(0,1);
            else
                Prob.Mode = 0;
                c_hat = nlp_c(x_hat, Prob, varargin{:});
            end
            z=[Ax;c_hat];
            cErr=max(z-bu,bl-z);
            if DEBUG
                xprinte(cErr,'cErr:')
            end
            hJ_hat = inf;
            % We should run a sequence of SOC steps. But how shall we update? The SOC
            % step is described on page 395 in the book "Practical Methods of
            % Optimization" by Fletcher.
            IterOne = 1;
            GOON=1;
            while GOON
                % Solve problem QP2 (page 7)
                ProbQP.x_0=zeros(n,1);

                Ad=ProbQP.A*d_hat;
                ProbQP.b_L=bl-z + Ad;
                ProbQP.b_U=bu-z + Ad;
                % optimization around x == 0, bounds from Prob.x_L,Prob.x_U
                ProbQP.x_L  = max(-min(rho)*ones(n,1),Prob.x_L-x_k);
                ProbQP.x_U  = min( min(rho)*ones(n,1),Prob.x_U-x_k);
                if DEBUG
                    fprintf('SOC Phase 2 %d\n',IterOne);
                    xprint(x_k,'x:');
                end


                if DEBUG
                    disp('qpPhase1: SEQUENCE OF SOC STEPS!!!')
                end
                ResultQP = tomRunMini(SolverQP,ProbQP);

                d_hat=ResultQP.x_k;
                l_hat=ResultQP.v_k;
                SOCerr=ResultQP.ExitFlag;

                if DEBUG
                    disp('qpPhase1: SOC2 result')
                    disp(SOCerr)
                    xprinte(Ad,'Ad:');
                    disp('bl z bu')
                    disp([bl z bu])
                    disp('x_L  x_U')
                    disp([ProbQP.x_L ProbQP.x_U])
                    disp('b_L  b_U')
                    disp([ProbQP.b_L ProbQP.b_U])
                    PrintMatrix(ProbQP.A,'A;dc:')
                    xprinte(d_hat,'d_hat:');
                    pause
                end

                % Criterion 2, page 8
                if SOCerr > 0
                    accept = 0;
                    r=NaN;

                    cErr=max(z-bu,bl-z);

                    if IterOne
                        x_hat = x_k;
                        dc_hat = dc_k;
                        % c_hat and c_exp computed
                        if sum(~J0) > 0
                            hJ_hat   = norm(max(0,cErr(~J0)),1);
                        else
                            hJ_hat   = 0;
                        end
                        if sum(J0) > 0
                            hT_hat = norm(max(0,cErr(J0)),1);
                        else
                            hT_hat = 0;
                        end
                        %hJ_hat = norm(max(0,cErr(J~=0)),1);
                        %hT_hat = norm(max(0,cErr(J==0)),1);
                    end
                    if DEBUG
                        disp('DID NOT CONVERGE WITH SOC STEP');
                        disp('TRY TO GO ON in qpPhase1');
                        pause
                    end
                    GOON=0;
                    %return                   % Terminate function
                else
                    IterOne = 0;

                    % Update parameters (x_hat = accepted point)
                    hJ_hatm1 = hJ_hat;
                    x_hat = x_k + d_hat;
                    % Adjust solution inside lower and upper variable bound
                    x_hat = max(Prob.x_L,min(Prob.x_U,x_hat));
                    if mA > 0, Ax=A*x_hat; end
                    if isempty(c)
                        c_hat=zeros(0,1);
                    else
                        Prob.Mode = 0;
                        c_hat = nlp_c(x_hat, Prob, varargin{:});
                    end
                    z=[Ax;c_hat];
                    cErr=max(z-bu,bl-z);
                    if DEBUG
                        xprinte(cErr,'cErr:')
                    end

                    if sum(~J0) > 0
                        hJ_hat   = norm(max(0,cErr(~J0)),1);
                    else
                        hJ_hat   = 0;
                    end
                    if sum(J0) > 0
                        hT_hat = norm(max(0,cErr(J0)),1);
                    else
                        hT_hat = 0;
                    end
                    %hJ_hat = norm(max(0,cErr(J~=0)),1);
                    %hT_hat = norm(max(0,cErr(J==0)),1);
                    r = hJ_hat/hJ_hatm1;          % SOC convergence rate

                    % Check if (hJ_hat,hT_hat) is acceptable to the filter
                    accept = CheckFilter1(FILTER, hJ_hat, hT_hat, cTol);

                    % Criterion 1, 3 and 4 on page 8
                    if DEBUG
                        disp(accept)
                        disp(r)
                        disp(hJ_hat)
                    end
                    if accept| r >= 0.25 | hJ_hat <= cTol
                        %return                   % Terminate function
                        GOON=0;
                    end
                end
            end

            % End of SOC2 part -------------------------------------

            if accept % (hJ_hat,hT_hat) is acceptable to the filter
                % Accept x_hat as new point
                d_k = x_hat - x_k;
                x_k = x_hat;
                % Do not update l_k (OK?)
                % Do not update delta_q and my (OK?)
                c_k = c_hat;
                % HKH ADDED
                if m==0
                    dc_k=zeros(0,1);
                else
                    Prob.Mode = 1;
                    dc_k = nlp_dc(x_k, Prob, varargin{:});
                end
                if DEBUG
                    disp('diff deriv')
                    disp(dc_k-dc_hat)
                end
                hJ_k = hJ_hat;
                hT_k = hT_hat;

                % Add (hJ_hat,hT_hat) to filter and remove any
                % points that are dominated by (hJ_hat,hT_hat)
                FILTER = ChangeFilter1(FILTER, hJ_k, hT_k, delta_q, my);

                % Possibly increase trust-region radius (item 2, page 16)
                if r <= 0.1
                    rho = 2*rho;
                end

            else
                if mA > 0, Ax=A*x_k; end
                z=[Ax;c_k];
                [h_k cErr] = ConstraintError(z-bu,bl-z,Alg);

                if Jsum~=0 & ( ~all(J==Jold) & h_k < u )
                    % Accept the best SOC step (set x_k=x_hat)
                    d_k = x_hat - x_k;
                    x_k = x_hat;
                    % Do not update l_k (OK?)
                    % Do not update delta_q and my (OK?)
                    c_k = c_hat;
                    dc_k = dc_hat;
                    hJ_k = hJ_hat;
                    hT_k = hT_hat;

                    % Remove all blocking entries that dominates
                    % (hJ_hat,hT_hat) from the filter (Section 3.4)
                    %FILTER = FILTER(find(FILTER(:,1) + xTol > hJ_k | ...
                    %         FILTER(:,2) + cTol > hT_k),:);
                    FILTER = FILTER(find(FILTER(:,1) > hJ_k | ...
                        FILTER(:,2) > hT_k),:);

                    % Add (hJ_hat,hT_hat) to filter? No, not if bellow works

                    % Reduce the upper bound u, add (-inf,u) to filter
                    % and and remove points dominated by (-inf,u)
                    u = max(h_k,u/10);
                    FILTER = ChangeFilter1(FILTER, -inf, u, 1, 1000000);

                else
                    % Reject the step (x_kp1=x_k, i.e. do not
                    % change x_k) Reduce the trust region radius
                    % (item 1, page 16)
                    if isempty(d_k)
                        rho = rho/2;
                        if DEBUG
                            disp('REDUCE rho, Reject SOC steps, keep x_k')
                            disp(rho)
                        end
                    else
                        % Safeguard update of rho
                        rho = min(rho,max(0.01*rho,norm(d_k,inf)))/2;
                        if DEBUG
                            disp('REDUCE rho, using previously accepted qpPhase1 step d_k ')
                            disp(d_k)
                            norm(d_k,inf)
                            disp(rho)
                        end
                    end
                end
            end
        end
    end
    Iter = Iter + 1;
    Jold = J;
    Jsum=sum(J);
end

ResultQ1.x_k = x_k;
ResultQ1.ExitFlag = 1;
if PriLev > 1
    fprintf('qpPhase1: Trust region rho too small. Failure!\n');
end
% ------------------------------
function Text = ExitText(ExitFlag,Inform)
% ------------------------------

switch  ExitFlag
    case 0
        switch  Inform
            case 1
                Text = 'Iteration points are close';
            case 2
                Text = 'Small search direction';
            case 3
                Text = str2mat('Function value below given estimate' ...
                    ,'Restart with lower Prob.f_Low if min not reached');
            case 4
                Text = 'Projected gradient small';
            case 10
                Text = 'Karush-Kuhn-Tucker conditions fulfilled';
        end
        Text = str2mat('Optimal point found',Text);
    case 1
        Text = 'Infeasible problem?';
    case 2
        Text = 'Maximal number of iterations reached';
    case 3
        Text = 'No progress in either f(x) or constraint reduction';
end


% MODIFICATION LOG:
%
% 981013  hkh  Safeguard code if c or dc empty. Also set empty as zeros(0,1)
% 981021  hkh  Check if qpopt is present, otherwise use qpSolve
% 981026  hkh  Use optParam.Method to make choice of QPOPT or qpSolve
%              Avoid infinite loop by setting ProbQP.QP.Phase1=0 in call to
%              qpSolve
% 981029  hkh  Trust region algorithm wrong, rho=0. Changed algorithm.
%              Corrected and safeguarded update of my.
%              Changed to avoid identity matrices in QP constraints
%              Avoid unnecessary loops building A_QP,b_QP,c_QP
% 981101  hkh  Big revision
% 981102  hkh  Use cTol in tests if constraints are violated
%              Remove eps setting, use cTol for test of hJ_hat
% Change test: if accept & ~all(d_kp1==0) avoiding d_kp1==0
% 990615  hkh  Change test on J changed, using Jold
% 001101  hkh  Use qld as default QP solver for v2.1 and v3.0 /MINI
% 020409  hkh  Add usage of Prob.nState and Prob.Mode
% 020820  ago  Change for Matlab 6.5 logical handling
% 030406  hkh  Lower-Upper bound bugs in qpPhase1 fixed.
% 030406  hkh  Add Lower-Upper bound adjustments for sharp simple bounds
% 031201  hkh  Revise algorithm selection for use of MAD and ADMAT
% 040111  hkh  Change call to inisolve, call NameAlg first, define mLin
% 040125  hkh  Set field mLin in ProbFP
% 041125  hkh  Norm was missing in KT=norm((g_k-v_k(1:n))/myMax);
%              Crucial - solver stopped too early
% 041126  hkh  Add test on proj.grad. not too high, if to accept small d_K
% 041126  hkh  Add test on |d_k| not too high, if to accept small proj.grad.
% 060814  med  FUNCS used for callbacks instead
% 060818  hkh  Use Prob.f_Low instead of optParam.eps_absf for convergence test
% 061212  med  ADMAT removed
% 070715  med  finite switched to isfinite
% 070907  hkh  SolverQP/FP picked from list, 1st with license, avoid GetSolver
% 080310  hkh  Must set ProbFP.QP.c  = zeros(n,1); due to system change in lp_f
% 080607  hkh  Use tomRun, not tomSolve
% 080607  med  Switched to tomRunMini
% 080911  med  Updated help to PriLevOpt