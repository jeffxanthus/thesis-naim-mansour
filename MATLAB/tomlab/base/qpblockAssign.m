% qpblockAssign implements the TOMLAB (TQ) format for nonlinearly
% constrained factorized QP problems.
%
% qpblockAssign sets the variables normally needed for an optimization in
% the TOMLAB structure Prob.
%
% It is specifically designed so that the user may give a sparse
% matrix and several additional blocks.
%
% -----------------------------------------------------
%
% Nonlinear constrained QP minimization problem:
%
%  Case 1:
%
%  The quadratic objective can be separated into main sparse matrix F and
%  one or more additional two part blocks.
%
%        min    0.5 * x' * F * x + d' * x +
%               0.5 * x' * Fb.out' * Fb.inn * Fb.out * x (for i=1...p)
%
%  Case 2:
%
%  The quadratic objective can be separated into a main sparse matrix F and
%  two outer blocks.
%
%        min    0.5 * x' * F * x + d' * x +
%               0.5 * x' * Fb.out' * Fb.out * x (for i=1...p)
%
%         x
%        s/t   x_L <=   x  <= x_U       x in R^n
%              b_L <= A x  <= b_U
%              c_L <= c(x) <= c_U
%
%  Case 3:
%
%  Case number 1 above, but F is supplied as a nx1 vector (diagonal)
%
%  Case 4:
%
%  Case number 2 above, but F is supplied as a nx1 vector (diagonal)
%
% Linear    equality equations: Set b_L==b_U
% Nonlinear equality equations: Set c_L==c_U
% Fixed     variables:          Set x_L==x_U. Both x_L and x_U must be finite.
%
% -----------------------------------------------------
%
% Syntax of qpblockAssign:
%
% function Prob = qpblockAssign(F, Fb, d, x_L, x_U, Name, x_0, ...
%                 A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U, ...
%                 fLowBnd, x_min, x_max, f_opt, x_opt);
%
% INPUT (Call with at least four parameters)
%
% F           The matrix F in 0.5 x' F x in the objective function
%
% Fb          A struct with the blocks for additional part of objective
%             function. For example:
%
%             Case 1:
%
%             Fb(1).out = ones(20,n)
%             Fb(1).inn = ones(20,20)
%
%             Fb(2).out = ones(30,n)
%             Fb(2).inn = ones(30,30)
%
%             Case 2:
%
%             Fb(1).out = ones(20,n)
%             Fb(2).out = ones(30,n)
%
%             Case 3:
%
%             Same as case 1, but F is a column vector (nx1).
%
%             Case 4:
%
%             Same as case 2, but F is a column vector (nx1).
%
% d           The vector d in d'x in the objective function
%
% x_L         Lower bounds on parameters x. If [] set as a nx1 -Inf vector.
% x_U         Upper bounds on parameters x. If [] set as a nx1  Inf vector.
%
% Note:       The number n of the unknown variables x are taken as
%             max(length(x_L),length(x_U),length(x_0))
%             You must specify at least one of these with correct length,
%             then the others are given default values
%
%             The following parameters are optional, and problem type dependent
%             Set empty to get default value
%
% Name        The name of the problem (string)
% x_0         Starting values, default nx1 zero vector
%
% L I N E A R   C O N S T R A I N T S
% A           mA x n matrix A, linear constraints b_L <= A*x <= b_U. Dense or sparse
% b_L         Lower bound vector in linear constraints b_L <= A*x <= b_U.
% b_U         Upper bound vector in linear constraints b_L <= A*x <= b_U.
%
% N O N L I N E A R   C O N S T R A I N T S
% c           Name of function that computes the mN nonlinear constraints
% dc          Name of function that computes the constraint Jacobian mN x n
% d2c         Name of function that computes the second part of the
%             Lagrangian function (only needed for some solvers)
%             See the help gateway routine nlp_d2c for an explanation of d2c
%
% ConsPattern mN x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the constraint Jacobian and ones indicate values that
%             might be non-zero. Used when estimating the Jacobian numerically.
%             Estimated before solve, if Prob.LargeScale==1, ConsPattern==[]
%
% c_L         Lower bound vector in nonlinear constraints c_L <= c(x) <= c_U.
% c_U         Upper bound vector in nonlinear constraints c_L <= c(x) <= c_U.
%
% A D D I T I O N A L   P A R A M E T E R S
% fLowBnd A lower bound on the function value at optimum. Default -1E300
%         A good estimate is not critical. Use [] if not known at all.
%
% x_min   Lower bounds on each x-variable, used for plotting
% x_max   Upper bounds on each x-variable, used for plotting
% f_opt   Optimal function value(s), if known (Stationary points)
% x_opt   The x-values corresponding to the given f_opt, if known.
%         If only one f_opt, give x_opt as a 1 by n vector
%         If several f_opt values, give x_opt as a length(f_opt) by n matrix
%         If adding one extra column n+1 in x_opt, 0 is min, 1 saddle, 2 is max.
%         x_opt and f_opt is used in printouts and plots.
%
% Set the variable as empty if this variable is not needed for the particular
% kind of problem you are solving

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.7.0$
% Written Aug 21, 2006.    Last modified Dec 13, 2006.

function Prob = qpblockAssign(F, Fb, d, x_L, x_U, Name, x_0, ...
    A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U, ...
    fLowBnd, x_min, x_max, f_opt, x_opt)

if nargin < 21
    x_opt=[];
    if nargin < 20
        f_opt=[];
        if nargin < 19
            x_max=[];
            if nargin < 18
                x_min=[];
                if nargin < 17
                    fLowBnd=[];
                    if nargin < 16
                        c_U=[];
                        if nargin < 15
                            c_L=[];
                            if nargin < 14
                                ConsPattern=[];
                                if nargin < 13
                                    d2c=[];
                                    if nargin < 12
                                        dc=[];
                                        if nargin < 11
                                            c=[];
                                            if nargin < 10
                                                b_U=[];
                                                if nargin < 9
                                                    b_L=[];
                                                    if nargin < 8
                                                        A=[];
                                                        if nargin < 7
                                                            x_0=[];
                                                            if nargin < 6
                                                                Name=[];
                                                                if nargin < 5
                                                                    error('qpblockAssign requires at least five parameters');
end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end

n = max([length(x_L),length(x_U),length(x_0)]);
if ~isempty(d)
    if ~isempty(F)
        if size(F,1) ~= length(d)
            error('Illegal length of d or number of rows in F.');
        end
    end
end
if ~isempty(A)
    if n ~= size(A,2)
        error('Illegal number of columns in A.');
    end
end

mQP = size(Fb,2);

qpcase = 1;
if isfield(Fb,'inn')
    if isempty(Fb(1).inn)
        qpcase = 2;
    end
else
    qpcase = 2;
end

if qpcase == 1
    for i=1:mQP
        if size(Fb(i).out,2) ~= n
            error('Illegal number of rows in Fb.out for number %d',i);
        end
        if size(Fb(i).out,1) ~= size(Fb(i).inn,1)
            error('Size of Fb.inn and Fb.out invalid for number %d',i);
        end
        if size(Fb(i).inn,1) ~= size(Fb(i).inn,2)
            error('Number of columns and rows don''t match for Fb.inn number %d',i);
        end
    end
else
    for i=1:mQP
        if size(Fb(i).out,2) ~= n
            error('Number of rows not correct for Fb.out number %d',i);
        end
    end
end

if size(F,1) ~= size(F,2)
    qpcase = qpcase+2;
    if length(F) ~= n
        error('Illegal length of F, case 3.');
    else
        F=F(:);
    end
end

Prob         = ProbDef(1);
Prob.QP.Fb   = Fb;
Prob.QP.mQP  = mQP;
Prob.QP.case = qpcase;
Prob = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A);

if isempty(d)
    d = [];
elseif length(d) ~= n
    error('Input d does not have correct length');
end

Prob.QP.F = F;
Prob.QP.c = double(full(d(:)));

Prob.P     = 1;
Prob.N     = n;
Prob.Name  = deblank(Name);
if isempty(x_min)
    Prob.x_min = -1*ones(n,1);
else
    Prob.x_min = x_min;
end
if isempty(x_max)
    Prob.x_max = 1*ones(n,1);
else
    Prob.x_max = x_max;
end
Prob.f_opt = f_opt;
Prob.x_opt = x_opt;

Prob.ConsPattern  = ConsPattern;

if ~isempty(fLowBnd)
    Prob.f_Low = max(Prob.f_Low,fLowBnd);
end

mN = max(length(c_L),length(c_U));
Prob.mNonLin = mN;

if mN > 0
    Prob.probType = checkType('con');
else
    Prob.probType = checkType('qp');
end

Prob.probFile=0;

if mN > 0
    if isempty(c_L)
        Prob.c_L=-Inf*ones(mN,1);
    else
        Prob.c_L=full(double(c_L(:)));
    end
    if isempty(c_U)
        Prob.c_U=Inf*ones(mN,1);
    else
        Prob.c_U=full(double(c_U(:)));
    end
    if any(Prob.c_L>Prob.c_U)
        error('c_L and c_U have crossover values');
    end
end
if mN == 0 & ~isempty(c)
    fprintf('WARNING in conAssign!!! ');
    fprintf('Constraint c is given. But no lower or upper bounds.\n');
end

if isempty(c) & (~isempty(c_L) | ~isempty(c_U))
    error('Constraint c not given, but lower and/or upper bounds.');
end

% Set Print Level to 0 as default
Prob.PriLevOpt=0;

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print
if isempty(MAX_x)
    MAX_x=20;
end
if isempty(MAX_c)
    MAX_c=20;
end
if isempty(MAX_r)
    MAX_r=30;
end

Prob = tomFiles(Prob, 'qpblock_f', 'qpblock_g', 'qpblock_H', funch2str(c), funch2str(dc), funch2str(d2c));

% MODIFICATION LOG
%
% 060821  med  Written, based on qpconAssign
% 060822  hkh  Safe guard F when vector to be a column vector
% 060822  hkh  Added safe guard on length of d
% 060824  med  Set Prob.QP.c to [] if empty
% 060824  med  Removed setting HessPattern to spones(F)
% 060824  med  Setting all vectors to full double
% 061213  med  Moved most input checks to checkAssign