% addPwLinFun - Adds piecewise linear function to a TOMLAB MIP problem.
%
% There are two ways to call addPwLinFun:
%
%  Syntax 1:
%    function Prob = addPwLinFun(Prob, 1, type, var, funVar, point, slope,
%                                a, fa)
%
%  Syntax 2:
%    function Prob = addPwLinFun(Prob, 2, type, var, funVar, firstSlope,
%                                point, value, lastSlope)
%
% This function will make one already existing variable of the problem to
% be constrained equal to a piecewise linear function of another already
% existing variable in the problem. The independent variable must be
% bounded in both directions.
%
% The variable constrained to be equal to a piecewise linear function can
% be used like any other variable; in constraints or the objective
% function.
%
% Depending on how many segments the function consists of, a number of new
% variables and constraints are added to the problem.
%
% Increasing the upper bound (x_U) or decreasing the lower bound (x_L) of
% the independent variable after calling this function will ruin the
% piecewise linear function.
%
% If the problem is to be solved by CPLEX, set type = 'cplex' to enhance
% performance. Otherwise, let type = 'mip'. NOTICE! You can not solve a
% problem with another solver than CPLEX if type = 'cplex'.
%
% INPUTS:
%
%   Prob           The problem to add the function to.
%   input          Flag indicating syntax used.
%   type           A string telling whether to construct a general MIP
%                  problem or to construct an MIP problem only solvable by
%                  CPLEX. Possible values: 'mip', 'cplex'
%   var            The number of the variable on which the piecewise linear
%                  function depends. Must exist in the problem already.
%   funVar         The number of the variable which will be equal to the
%                  piecewise linear function. Must exist in the problem
%                  already.
%   firstSlope     Syntax 2 only. The slope of the piecewise linear
%                  function left of the first point, point(1).
%   point          An array of break points. Must be sorted. If two values
%                  occur twice, there is a step at that point. Length r.
%   slope          Syntax 1 only. An array of the slopes of the segments.
%                     slope(i) is the slope between point(i-1) and point(i).
%                     slope(1) is the slope of the function left of point(1).
%                     slope(r+1) is the slope of the function right of
%                      point(r).
%                  If points(i-1) == points(i), slope(i) is the height of
%                  the step.
%   value          Syntax 2 only. The values of the piecewise linear
%                  function at the points given in point. f(point(i)) =
%                  value(i). If point(i-1) == point(i), value(i-1) is the
%                  right limit of the value at the point, and value(i) is
%                  the left limit of the value at the point.
%   lastSlope      Syntax 2 only. The slope of the piecewise linear
%                  function right of the last point, point(r).
%   a, fa          Syntax 1 only. The value of the piecewise linear
%                  function at point a is equal to fa. f(a) = fa, that is.
%
% OUTPUT:
%
%   Prob           The new problem structure with the piecewise linear
%                  function added. New variables and linear constraints
%                  added. (MIP problem)

% Fredrik Hellman, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2007-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Jul 23, 2007.   Last modified Aug 13, 2009.

function Prob = addPwLinFun(Prob, input, type, var, funVar, varargin)

if input == 1 & nargin == 9
    point = varargin{1};
    slope = varargin{2};
    a     = varargin{3};
    fa    = varargin{4};

    % Check sizes. slope should always be one element larger than point
    if(length(point) ~= length(slope)-1)
        error('length(point) + 1 should be equal to length(slope)');
    end

    r = length(point);
    value = zeros(r,1);

    for i = 2:r
        if(point(i-1) == point(i))
            value(i) = value(i-1) + slope(i);
            if(point(i) == a)
                error('a must not be a step point.');
            end
        else
            value(i) = value(i-1) + (point(i)-point(i-1))*slope(i);
        end
        if point(i-1) <= a & point(i) > a
            lambda = (a-point(i-1))/(point(i)-point(i-1));
            tfa = (1-lambda)*value(i-1)+lambda*value(i);
        end
    end

    firstSlope = slope(1);
    lastSlope  = slope(end);

    if(a < point(1))
        tfa = value(1)+firstSlope*(a-point(1));
    elseif(a >= point(end))
        tfa = value(r)+lastSlope*(a-point(r));
    end

    % tfa is the value of f(a) with the current point/value.
    % It should be fa. We need to add fa-tfa to all value elements.

    value = value + fa - tfa;

    %point       point1    point2   point3   point4
    %slope slopeLeft, slope12, slope23, slope34, slopeRight
    %value       value1    value2   value3   value4
elseif(nargin == 9)
    firstSlope = varargin{1};
    point      = varargin{2};
    value      = varargin{3};
    lastSlope  = varargin{4};

    if(length(point) ~= length(value))
        error('length(point) should be equal to length(value)');
    end
else
    error('need exactly 9 arguments.');
end

point = point(:)';
value = value(:)';

if(any(diff(find(diff(point) == 0)) == 1))
    error('a maximum of two points can be equal');
end

if(~issorted(point))
    error('point must be sorted');
end

if(strcmp(type, 'cplex'))
    cplex = 1;
else
    cplex = 0;
end

if funVar < 1 | funVar > Prob.N
    error('funVar must be one of 1, ..., Prob.N');
end

if var < 1 | var > Prob.N
    error('var must be one of 1, ..., Prob.N');
end

if funVar == var
    error('funVar is equal to var');
end

N = Prob.N;

var_L = Prob.x_L(var);
var_U = Prob.x_U(var);

% Require bounds on the variable
if(~isfinite(var_U-var_L))
    error('The independent variable must be bounded.');
end

% Final variable order:
%
% [original variables, y_left, y_right, (lambda_1), lambda2, ..., lambda_r-1,
% (lambda_r), x_left, x_right, x_1, ..., x_s, p_left, p_right]

% Figure out where there are steps.

r       = length(point);
dpoint  = diff(point);
stepidx = find(dpoint==0);
if(cplex)
    nsets = length(stepidx)+1;
    s = nsets;
else
    nsegments = nnz(dpoint);
    s = r-1;
end

% The problem is formulated according to the formulation at page 11 in
% Integer and Combinatorial Optimization by Nemhause and Wolsey, 1988, with
% two modifications:
%
%  1. Special treatment of the first and last segment.
%  2. Support for discontinuous functions. (zero-width segments).

% Lambda constraints

if(cplex)
    % In case of cplex, we have only one constraint per set. The sum of
    % the lambdas belonging to a set is equal to the binary variable x
    % corresponding to the set.
    b_Llambda = zeros(nsets, 1);
    b_Ulambda = zeros(nsets, 1);
    Alambda   = sparse(nsets, N+2+r+2+s+2);
    setedges = [0; stepidx(:); r];
    for i = 1:nsets
        lambda_idx = N+2+((setedges(i)+1):setedges(i+1));
        x_idx      = N+2+r+2+i;
        Alambda(i,lambda_idx) = -1;
        Alambda(i,x_idx)      = 1;
    end
    nlambdacon = nsets;
else
    b_Llambda = zeros(r,1);
    b_Ulambda = Inf*ones(r,1);
    Alambda = [sparse(r, N+2), -speye(r,r), sparse(r, 2), [speye(s,s); sparse(1,s)] + [sparse(1,s); speye(s,s)], sparse(r, 2)];
    nlambdacon = r;
end
% sum(x_i) = 1, i = 0, ..., s+1

b_Lsegsum = 1;
b_Usegsum = 1;
Asegsum = [sparse(1, N+2+r), spones(ones(1, 2+(s))), sparse(1,2)];

% sum(x_i) = sum(lambda_j), i = 1, ..., s, j = 1, ..., r

b_Llssum = 0;
b_Ulssum = 0;
Alssum  = sparse([sparse(1, N+2), spones(ones(1, r)), sparse(1, 2), -spones(ones(1,(s))), sparse(1, 2)]);

% 0 = var - sum(lambda_i * a_i) -
%     + p_left - a_1 * x_left - p_right - a_r * x_right

b_Lvar  = 0;
b_Uvar  = 0;
Avar    = sparse([sparse(1, N+2), -point, -point(1), -point(end), sparse(1, (s)), 1, -1]);
Avar(var) = 1;

% Now the piecewise linear function itself
%
% f(y) = sum(lambda_i * a_i) - firstSlope * p_left + f(a_1) * x_left -
%        + lastSlope * p_right + f(a_r) * x_right
%
% f(y) is funVar, so we will add another constraint:
%
% 0 = sum(lambda_i * a_i) - firstSlope * p_left + f(a_1) * x_left -
%        + lastSlope * p_right + f(a_r) * x_right - f(y)

b_Lfunc  = 0;
b_Ufunc  = 0;
Afunc    = sparse([sparse(1,N+2), value, value(1), value(end), sparse(1, (s)), -firstSlope, lastSlope]);
Afunc(funVar) = -1;

b_Lnew  = [b_Llambda; b_Lsegsum; b_Llssum; b_Lvar; b_Lfunc];
b_Unew  = [b_Ulambda; b_Usegsum; b_Ulssum; b_Uvar; b_Ufunc];
Anew    = sparse([Alambda; Asegsum; Alssum; Avar; Afunc]);

% Now, remove some variables:
% In case of cplex, we need to remove the sets with only one lambda
% connected to it. These are zero-width sets.
%
% Othwerise, the segments that are not used, those that are in a step
% should be removed. The indices of these are givne in stepidx
% At the same time, update s = nsegments;

if(cplex)
    remlambdaidx = [];
    dsetedges = diff(setedges);
    zerosetidx = find(dsetedges == 1);
    % This can only occur in the ends.
    lambdafirst = 1;
    lambdalast  = setedges(end);
    if(any(zerosetidx == 1))
        remlambdaidx(end+1) = lambdafirst;
        setedges(1) = [];
        setedges = setedges-1;
    end
    if(any(zerosetidx == nsets))
        remlambdaidx(end+1) = lambdalast;
        setedges(end) = [];
    end
    Anew(:,N+2+r+2+zerosetidx) = [];
    nsets = nsets - length(zerosetidx);
    s = nsets;
else
    Anew(:,N+2+r+2+stepidx) = [];
    s = nsegments;
end

% The first and/or last lambda variables might be fixed. Remove them in that
% case.

% The rows of Alambda (the first part of Anew) that have a row sum that is
% equal to -1, corresponds to fixed fixed lambdas for both cplex and
% non-cplex. Remove those rows, and the corresponding variables.

if(~cplex)
    remlambdaidx = find(sum(Anew(1:nlambdacon,:),2) == -1);
end
nlambdacon = nlambdacon - length(remlambdaidx);

% Remove variables:
Anew(:,N+2+remlambdaidx) = [];

% There might be empty rows now (sum(x_i) = sum(lambda_i)) if there are
% no lambda_i or x_i.
emptyrowidx = find(sum(abs(Anew), 2) == 0);

% Remove constraints:
Anew(emptyrowidx,:)   = [];
b_Lnew(emptyrowidx)   = [];
b_Unew(emptyrowidx)   = [];

% We might now have fewer lambdas
nlambda = r-length(remlambdaidx);

% Note that y_right/left and p_right/left are shifted and defined in the
% interval [0, inf]. This is needed to make the bincont2lin not complain.
% The left-variables (y_left and p_left) are negative of their "actual"
% value as well.

% Then we have the constraints:
%  p_left  = x_left   * y_left
%  p_right = x_right  * y_right
% These are set through the bincont2lin utility at the end.
% Box bounds:
%           0 <= y_left   <= max(a_1-var_L, 0)
%           0 <= y_right  <= max(var_U-a_r, 0)
%           0 <= lambda_i <= inf  (i = 1,...,nlambda)
%           0 <= x_left   <= 1
%           0 <= x_right  <= 1
%           0 <= x_i      <= 1    (i = 1,...,s)
%           0 <= p_left   <= max(a_1-var_L, 0)
%           0 <= p_right  <= max(var_U-a_r, 0)

x_Lnew = [0, 0, zeros(1,nlambda), 0, 0, zeros(1, s), 0, 0]';
x_Unew = [max(point(1)-var_L, 0), max(var_U-point(end), 0), inf*ones(1,nlambda), 1, 1, ones(1, s), max(point(1)-var_L, 0), max(var_U-point(end), 0)]';

Nnew   = length(x_Lnew);

if(cplex)
    % Formulate the SOS2
    sos2 = struct('var', [], 'row', []);
    sos2(nsets).var = [];

    for i = 1:nsets
        sos2(i).var = N + 2 + ((setedges(i)+1):setedges(i+1));
        % The var row contains the point values as coefficients for the lambda
        % variables. The point values are sorted in the order of which the
        % variables in the SOS2 set should be ordered.
        sos2(i).row = Prob.mLin + nlambdacon + 3;
    end
end

% Update the Prob structure

Prob.N = Prob.N + Nnew;

if(~isempty(Prob.QP.c))
    Prob.QP.c = [Prob.QP.c(:); zeros(Nnew, 1)];
end

if(isfield(Prob.QP, 'F') && ~isempty(Prob.QP.F))
    Prob.QP.F = blkdiag(Prob.QP.F, zeros(Nnew));
end

if(isfield(Prob.QP, 'B') && ~isempty(Prob.QP.B))
    Prob.QP.B = [Prob.QP.B(:); zeros(Nnew, 1)];
end

Prob.x_L = [Prob.x_L(:); x_Lnew];
Prob.x_U = [Prob.x_U(:); x_Unew];

if(~isempty(Prob.A))
    Prob.A = [Prob.A, sparse(Prob.mLin, Nnew); Anew];
else
    Prob.A = Anew;
end
Prob.b_L = [Prob.b_L(:); b_Lnew];
Prob.b_U = [Prob.b_U(:); b_Unew];

Prob.mLin = size(Prob.A,1);

IntVarIdx = N+2+nlambda+(1:s+2)';
if(isempty(Prob.MIP.IntVars))
    % No intvar array. Create an index-array.
    Prob.MIP.IntVars = IntVarIdx;
elseif(length(Prob.MIP.IntVars) == N)
    % Logical array
    Prob.MIP.IntVars = [Prob.MIP.IntVars(:); zeros(Nnew,1)];
    Prob.MIP.IntVars(IntVarIdx) = 1;
else
    % An already existing index-array
    Prob.MIP.IntVars = [Prob.MIP.IntVars(:); IntVarIdx];
end

if(isempty(Prob.MIP.VarWeight))
elseif(length(Prob.MIP.VarWeight) == N)
    Prob.MIP.VarWeight = [Prob.MIP.VarWeight(:); zeros(Nnew,1)];
else
    error('Prob.MIP.VarWeight has invalid length');
end

if(cplex)
    if(~isempty(Prob.MIP.sos2))
        oldsos2 = Prob.MIP.sos2;
    else
        oldsos2 = [];
    end
    Prob.MIP.sos2 = [oldsos2(:)', sos2];
end

% Done updating
% Now add the bincont2lin constraints
Prob = bincont2lin(Prob, N+2+nlambda+2+s+[1;2], N+2+nlambda+[1;2], N+[1;2]);

% MODIFICATION LOG:
%
% 070723 frhe  Written
% 070725 frhe  First version complete.
% 070726 frhe  Piecewise linear function now equal to a variable
% 070726 frhe  Now use SOS2 sets for CPLEX
% 070727 frhe  Support for syntax 1 added
% 070806 med   zeros --> sparse for most locations
% 080903 frhe  Corrected minor bug
% 090813 med   mlint check