% function exitFlag = checkFuncs(Prob, Solver, PriLev)
%
% TOMLAB routine for verifying user supplied routines
%
% INPUT PARAMETERS
% Prob         Problem structure created with assign routine.
% Solver       Solver that will be used. For example 'knitro' (default).
% PriLev       0 - suppress warnings (info), 1 - full printing (default)
%
% OUTPUT PARAMETERS
% exitFlag     0 if no errors.

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Feb 13, 2005.   Last modified Aug 14, 2006.

function exitFlag = checkFuncs(Prob, Solver, PriLev)

if nargin < 1
    error('checkFuncs requires at least one parameter, Prob');
end

if nargin < 3
    PriLev = 1;
    if nargin < 2
        Solver='knitro';
    end
end

% Run ProbCheck to fix problems.

Prob = ProbCheck(Prob, Solver);
FUNCS = Prob.FUNCS;
x_0 = Prob.x_0;
x_L = Prob.x_L;
x_U = Prob.x_U;
N   = Prob.N;

c_L = Prob.c_L;
c_U = Prob.c_U;

A = Prob.A;
b_L = Prob.b_L;
b_U = Prob.b_U;

% CHECK BOUNDS, x_0

if isempty(x_0)
    x_0 = zeros(N,1);
elseif length(x_0) ~= N
    fprintf('Length of x_0 %d, should be %d\n',length(x_0),N);
    fprintf('Illegal length of x_0');
    error('Cannot proceed with more testing until fixed.');
else
    x_0 = x_0(:);
end

if isempty(x_L)
    x_L = -Inf*ones(N,1);
elseif length(x_L) ~= N
    fprintf('Length of x_L %d, should be %d\n',length(x_L),N);
    error('Illegal length of x_L');
else
    x_L = x_L(:);
end

if isempty(x_U)
    x_U = Inf*ones(N,1);
elseif length(x_U) ~= N
    fprintf('Length of x_U %d, should be %d\n',length(x_U),N);
    error('Illegal length of x_U');
else
    x_U = x_U(:);
end

% CHECK ON LINEAR CONSTRAINTS

[mA,mN] = size(A);

if mA > 0
    Prob.A = A;
    % MAKE SURE A CORRECT.
    if mN ~= N
        fprintf('Number of variables %d\n',N);
        fprintf('Number of columns in linear constraint matrix A %d\n',mN);
        fprintf('These lengths should be the same, check input.\n');
        fprintf('Illegal number of columns in A.\n \n');
    end
    % MAKE SURE b_L CORRECT
    if size(b_L,1) ~= mA
        fprintf('Number of linear constraints %d\n',mA);
        fprintf('Number of rows in lower bounds on linear constraint %d\n',size(b_L,1));
        fprintf('These lengths should be the same, check input.\n');
        fprintf('Illegal number of rows in b_L.\n \n');
    end
    % MAKE SURE b_U CORRECT
    if size(b_U,1) ~= mA
        fprintf('Number of linear constraints %d\n',mA);
        fprintf('Number of rows in upper bounds on linear constraint %d\n',size(b_U,1));
        fprintf('These lengths should be the same, check input.\n');
        fprintf('Illegal number of rows in b_U.\n \n');
    end
    % MAKE SURE NO DOUBLE SIDED CONSTRAINTS DOUBLE SIDED. FIND ROW THAT ARE
    % SAME
    for i=1:mA
        for j=i+1:mA
            if all(A(i,:) == A(j,:))
                fprintf('There are duplicate constraints in A on rows %d\n',i);
                fprintf('and %d. Remove extra constraints or move bounds\n',j);
                fprintf('in b_L and b_U.\n');
                fprintf('The format for TOMLAB is double-sided, b_L <= Ax <= b_U.\n');
            end
        end
    end
end

% Start checking which files exist. Problem types are LS (0), (MI)NLP (1)
% and SIM (2)

% START WITH LS.

if ~isempty(FUNCS.r)
    probType = 0;
elseif ~isempty(FUNCS.fc)
    probType = 2;
else
    probType = 1;
end

if ~isempty(x_0)
    x_0  = max(x_L(1:N),min(x_0,x_U(1:N)));
else
    if PriLev == 1
        fprintf('Prob.x_0 is empty, setting to (Prob.x_L+Prob.x_U)/2\n');
    end
    x_L = max(x_L,ones(length(x_L),1)*(-1e6));
    x_U = min(x_U,ones(length(x_L),1)*(+1e6));
    x_0 = (x_L+x_U)/2;
end

if probType == 0
    % CHECKS FOR LS PROBLEMS
    narginr  = xnargin(FUNCS.r);
    nargoutr = xnargout(FUNCS.r);

    % MAKE SURE CORRECT NUMBER OF INPUTS TO FUNCS.r
    if ~any(narginr == [1 2])
        fprintf('Incorrect number of inputs to user routine %s. \n',FUNCS.r);
        fprintf('The inputs should be x (and optionally Prob). \n');
        error('Cannot proceed with more testing until fixed.');
    end

    % PERFORM SIZE CHECK ON RESIDUAL ROUTINE

    if narginr > 1
        resid = feval(FUNCS.r, x_0, Prob);
    else
        resid = feval(FUNCS.r, x_0);
    end

    rlen  = length(resid);

    if isempty(Prob.LS.y)
        fprintf('y, the fit vector is empty. This is needed for memory allocation. \n');
        error('The vector is not used for nonlinear LS problems, set to zeros().');
    end

    ylen = length(Prob.LS.y);

    if ylen ~= rlen
        fprintf('Incorrect length of residuals in %s \n',FUNCS.r);
        fprintf('or of fit vector y. \n');
        error('Cannot proceed with more testing until fixed.');
    end

    if size(Prob.LS.y,2) > 1
        fprintf('The fit vector y is incorrect, it should be \n');
        fprintf('a column vector. \n');
        error('Cannot proceed with more testing until fixed.');
    end

    if size(resid,2) > 1
        fprintf('The residual vector in user routine %s is incorrect, \n', FUNCS.r);
        fprintf('it should be a column vector, allocate with r = zeros(m,1). \n');
        error('Cannot proceed with more testing until fixed.');
    end

    % MAKE SURE CORRECT NUMBER OF OUTPUTS FROM FUNCS.r
    if nargoutr ~= 1
        fprintf('Incorrect number of outputs from user routine %s. \n',FUNCS.r);
        fprintf('The output should be a residual vector. \n');
        error('Cannot proceed with more testing until fixed.');
    end

    if ~isempty(FUNCS.J)
        narginJ  = xnargin(FUNCS.J);
        nargoutJ = xnargout(FUNCS.J);

        % MAKE SURE CORRECT NUMBER OF INPUTS TO FUNCS.J
        if narginJ ~= narginr
            fprintf('The number of inputs to user routines %s and \n',FUNCS.r);
            fprintf('%s do not match. \n',FUNCS.J);
            error('Cannot proceed with more testing until fixed.');
        end

        % MAKE SURE CORRECT NUMBER OF OUTPUTS FROM FUNCS.J
        if nargoutJ ~= 1
            fprintf('Incorrect number of outputs from user routine %s. \n',FUNCS.J);
            fprintf('The output should be a Jacobian matrix. \n');
            error('Cannot proceed with more testing until fixed.');
        end

        % PERFORM SIZE CHECK ON JACOBIAN ROUTINE

        if narginJ > 1
            jac = feval(FUNCS.J, x_0, Prob);
        else
            jac = feval(FUNCS.J, x_0);
        end

        [mJac,nJac] = size(jac);

        if mJac == N & nJac == rlen & PriLev == 1
            fprintf('Warning - The Jacobian returned from user routine %s. \n',FUNCS.J);
            fprintf('may be transposed. The size should be mxn. \n');
            fprintf('J = J'' will fix the problem. \n');
        end

        if mJac ~= rlen
            fprintf('Incorrect number of rows from user routine %s. \n',FUNCS.J);
            fprintf('There are %d rows, but should be %d. \n',mJac,rlen);
            error('Cannot proceed with more testing until fixed.');
        end

        if nJac ~= N
            fprintf('Incorrect number of columns from user routine %s. \n',FUNCS.J);
            fprintf('There are %d columns, but should be %d. \n',nJac,N);
            error('Cannot proceed with more testing until fixed.');
        end
    end

    if ~isempty(FUNCS.d2r)
        nargind2r  = xnargin(FUNCS.d2r);
        nargoutd2r = xnargout(FUNCS.d2r);

        % MAKE SURE CORRECT NUMBER OF INPUTS TO FUNCS.d2r
        if nargind2r ~= narginr + 2
            fprintf('The number of inputs to user routines %s and \n',FUNCS.d2r);
            fprintf('%s do not match. \n',FUNCS.d2r);
            error('Cannot proceed with more testing until fixed.');
        end

        % MAKE SURE CORRECT NUMBER OF OUTPUTS FROM FUNCS.d2r
        if nargoutd2r ~= 1
            fprintf('Incorrect number of outputs from user routine %s. \n',FUNCS.d2r);
            fprintf('The output should be a Jacobian matrix. \n');
            error('Cannot proceed with more testing until fixed.');
        end

        % PERFORM SIZE CHECK ON D2R ROUTINE - NOT DONE.
    end
end

if probType == 1
    % CHECKS FOR (MI)NLP PROBLEMS
    narginf  = xnargin(FUNCS.f);
    nargoutf = xnargout(FUNCS.f);

    % MAKE SURE CORRECT NUMBER OF INPUTS TO FUNCS.f
    if ~any(narginf == [1 2])
        fprintf('Incorrect number of inputs to user routine %s. \n',FUNCS.f);
        fprintf('The inputs should be x (and optionally Prob). \n');
        error('Cannot proceed with more testing until fixed.');
    end

    % MAKE SURE CORRECT NUMBER OF OUTPUTS FROM FUNCS.f
    if nargoutf ~= 1
        fprintf('Incorrect number of outputs from user routine %s. \n',FUNCS.f);
        fprintf('The output should be a scalar objective. \n');
        error('Cannot proceed with more testing until fixed.');
    end

    % PERFORM SIZE CHECK ON FUNCTION ROUTINE

    if narginf > 1
        funcval = feval(FUNCS.f, x_0, Prob);
    else
        funcval = feval(FUNCS.f, x_0);
    end

    flen  = length(funcval);

    if length(flen) ~= 1
        fprintf('The objective function in user routine %s is incorrect, \n', FUNCS.f);
        fprintf('it should be a scalar (real values only). \n');
        error('Cannot proceed with more testing until fixed.');
    end

    if ~isempty(FUNCS.g)
        narging  = xnargin(FUNCS.g);
        nargoutg = xnargout(FUNCS.g);

        % MAKE SURE CORRECT NUMBER OF INPUTS TO FUNCS.g
        if narging ~= narginf
            fprintf('The number of inputs to user routines %s and \n',FUNCS.g);
            fprintf('%s do not match. \n',FUNCS.f);
            error('Cannot proceed with more testing until fixed.');
        end

        % MAKE SURE CORRECT NUMBER OF OUTPUTS FROM FUNCS.g
        if nargoutg ~= 1
            fprintf('Incorrect number of outputs from user routine %s. \n',FUNCS.g);
            fprintf('The output should be a gradient vector. \n');
            error('Cannot proceed with more testing until fixed.');
        end

        % PERFORM SIZE CHECK ON GRADIENT ROUTINE

        if narging > 1
            gvec = feval(FUNCS.g, x_0, Prob);
        else
            gvec = feval(FUNCS.g, x_0);
        end

        [grows,gcols] = size(gvec);

        if gcols ~= 1
            fprintf('Incorrect number of columns from user routine %s. \n',FUNCS.g);
            fprintf('There should only be 1 column with %d rows. \n', N);
            error('Cannot proceed with more testing until fixed.');
        end

        if grows ~= N
            fprintf('Incorrect number of rows from user routine %s. \n',FUNCS.g);
            fprintf('There should only be %d rows. \n', N);
            error('Cannot proceed with more testing until fixed.');
        end
    end

    if ~isempty(FUNCS.H)
        narginH  = xnargin(FUNCS.H);
        nargoutH = xnargout(FUNCS.H);

        % MAKE SURE CORRECT NUMBER OF INPUTS TO FUNCS.H
        if narginH ~= narginf
            fprintf('The number of inputs to user routines %s and \n',FUNCS.H);
            fprintf('%s do not match. \n',FUNCS.f);
            error('Cannot proceed with more testing until fixed.');
        end

        % MAKE SURE CORRECT NUMBER OF OUTPUTS FROM FUNCS.H
        if nargoutH ~= 1
            fprintf('Incorrect number of outputs from user routine %s. \n',FUNCS.H);
            fprintf('The output should be a Hessian matrix. \n');
            error('Cannot proceed with more testing until fixed.');
        end

        % PERFORM SIZE CHECK ON HESSIAN ROUTINE

        if narginH > 1
            Hmat = feval(FUNCS.H, x_0, Prob);
        else
            Hmat = feval(FUNCS.H, x_0);
        end

        [Hrows,Hcols] = size(Hmat);

        if Hrows ~= Hcols
            fprintf('Incorrect Hessian from user routine %s. \n',FUNCS.H);
            fprintf('The rows and columns are not the same. \n');
            fprintf('It should be a %d by %d matrix. \n', N, N);
            error('Cannot proceed with more testing until fixed.');
        end

        if Hrows ~= N
            fprintf('Incorrect Hessian from user routine %s. \n',FUNCS.H);
            fprintf('It should be a %d by %d matrix. \n', N, N);
            error('Cannot proceed with more testing until fixed.');
        end

        % MAKE SURE HESSIAN IS PERFECTLY SYMMETRIC.

        if any(find(Hmat-Hmat'))
            fprintf('The Hessian of the objective is not symmetric \n')
            fprintf('You must either correct your computations so that \n')
            fprintf('the difference H-H'' is of the order of 1E-16 at most \n')
            fprintf('or make the Hessian symmetric by doing the following trick: \n')
            fprintf('H = 0.5*(H+H''); \n')
            error('Cannot proceed with more testing until fixed.');
        end
    end
end

if isempty(FUNCS.c) & (length(c_L) > 0 | length(c_U) > 0)
    fprintf('There are bounds in either Prob.c_L or Prob.c_U, \n');
    fprintf('but no user routine for the constraints has been specified. \n');
    error('Cannot proceed with more testing until fixed.');
end

if ((probType == 0 | probType == 1) & ~isempty(FUNCS.c)) & (length(c_L) > 0 | length(c_U) > 0)

    if probType == 0
        nargincorr  = xnargin(FUNCS.r);
    else
        nargincorr  = xnargin(FUNCS.f);
    end

    narginc  = xnargin(FUNCS.c);
    nargoutc = xnargout(FUNCS.c);

    % MAKE SURE CORRECT NUMBER OF INPUTS TO FUNCS.c
    if narginc ~= nargincorr
        fprintf('Incorrect number of inputs to user routine %s. \n',FUNCS.c);
        fprintf('The inputs should be x (and optionally Prob). \n');
        error('Cannot proceed with more testing until fixed.');
    end

    % MAKE SURE CORRECT NUMBER OF OUTPUTS FROM FUNCS.c
    if nargoutc ~= 1
        fprintf('Incorrect number of outputs from user routine %s. \n',FUNCS.c);
        fprintf('There shoudl only be one output. Nonlinear constraints  \n');
        fprintf('in TOMLAB are of the form c_L < c(x) < c_U \n');
        fprintf('Equality constraints are introduced by setting c_L(i) == c_U(i) \n');
        error('Cannot proceed with more testing until fixed.');
    end

    % PERFORM SIZE CHECK ON CONSTRAINT ROUTINE

    if narginc > 1
        convec = feval(FUNCS.c, x_0, Prob);
    else
        convec = feval(FUNCS.c, x_0);
    end

    clen  = length(convec);

    if clen ~= max(length(c_L),length(c_U))
        fprintf('The consstraint vector returned by user routine %s. \n', FUNCS.c);
        fprintf('is of incorrect length. If Prob.c_L and Prob.c_U are correct \n');
        fprintf('the length should be %d. \n',max(length(c_L),length(c_U)));
        error('Cannot proceed with more testing until fixed.');
    end

    if size(convec,2) > 1
        fprintf('The constraint vector in user routine %s is incorrect, \n', FUNCS.c);
        fprintf('it should be a column vector, allocate with c = zeros(m,1). \n');
        error('Cannot proceed with more testing until fixed.');
    end

    if ~isempty(FUNCS.dc)
        nargindc  = xnargin(FUNCS.dc);
        nargoutdc = xnargout(FUNCS.dc);

        % MAKE SURE CORRECT NUMBER OF INPUTS TO FUNCS.dc
        if nargindc ~= nargincorr
            fprintf('The number of inputs to user routine %s is \n',FUNCS.dc);
            fprintf('incorrect. There should be %d inputs \n',nargincorr);
            fprintf('They should be x (and optionally Prob). \n');
            fprintf('Prob has to be given if given for other user routines. \n');
            fprintf('This includes TOMLAB interface routines for which Prob is given. \n');
            error('Cannot proceed with more testing until fixed.');
        end

        % MAKE SURE CORRECT NUMBER OF OUTPUTS FROM FUNCS.dc
        if nargoutdc ~= 1
            fprintf('Incorrect number of outputs from user routine %s. \n',FUNCS.dc);
            fprintf('The output should be a constraint Jacobian matrix. \n');
            error('Cannot proceed with more testing until fixed.');
        end

        % PERFORM SIZE CHECK ON JACOBIAN ROUTINE

        if nargindc > 1
            cjac = feval(FUNCS.dc, x_0, Prob);
        else
            cjac = feval(FUNCS.dc, x_0);
        end

        [mcJac,ncJac] = size(cjac);

        if mcJac == N & ncJac == clen & PriLev == 1
            fprintf('Warning - The constraint Jacobian returned from user routine %s. \n',FUNCS.dc);
            fprintf('may be transposed. The size should be mxn. \n');
            fprintf('dc = dc'' will fix the problem. \n');
        end

        if mcJac ~= clen
            fprintf('Incorrect number of rows from user routine %s. \n',FUNCS.dc);
            fprintf('There are %d rows, but should be %d. \n',mcJac,clen);
            error('Cannot proceed with more testing until fixed.');
        end

        if ncJac ~= N
            fprintf('Incorrect number of columns from user routine %s. \n',FUNCS.dc);
            fprintf('There are %d columns, but should be %d. \n',ncJac,N);
            error('Cannot proceed with more testing until fixed.');
        end
    end

    if ~isempty(FUNCS.d2c)
        nargind2c  = xnargin(FUNCS.d2c);
        nargoutd2c = xnargout(FUNCS.d2c);

        % MAKE SURE CORRECT NUMBER OF INPUTS TO FUNCS.d2c
        if nargind2c ~= nargincorr+1
            fprintf('The number of inputs to user routine %s is \n',FUNCS.d2c);
            fprintf('incorrect. There should be %d inputs \n',nargincorr+1);
            fprintf('They should be x (and optionally Prob). \n');
            fprintf('Prob has to be given if given for other user routines. \n');
            fprintf('This includes TOMLAB interface routines for which Prob is given. \n');
            error('Cannot proceed with more testing until fixed.');
        end

        % MAKE SURE CORRECT NUMBER OF OUTPUTS FROM FUNCS.d2r
        if nargoutd2c ~= 1
            fprintf('Incorrect number of outputs from user routine %s. \n',FUNCS.d2c);
            fprintf('The output should be d2c = sum_i=1:m  lam(i) * d2c(i) / dx^2. \n');
            error('Cannot proceed with more testing until fixed.');
        end

        % PERFORM SIZE CHECK ON D2C ROUTINE.

        if nargind2c > 2
            d2chess = feval(FUNCS.d2c, x_0, ones(clen,1), Prob);
        else
            d2chess = feval(FUNCS.dc, x_0, ones(clen,1));
        end

        [d2crows,d2ccols] = size(d2chess);

        if d2crows ~= d2ccols
            fprintf('Incorrect d2c from user routine %s. \n',FUNCS.d2c);
            fprintf('The rows and columns are not the same. \n');
            fprintf('It should be a %d by %d matrix. \n', N, N);
            error('Cannot proceed with more testing until fixed.');
        end

        if d2crows ~= N
            fprintf('Incorrect d2c from user routine %s. \n',FUNCS.d2c);
            fprintf('It should be a %d by %d matrix. \n', N, N);
            error('Cannot proceed with more testing until fixed.');
        end

        % MAKE SURE D2C IS PERFECTLY SYMMETRIC.

        if any(find(d2chess-d2chess'))
            fprintf('The Hessian of the constraints is not symmetric \n')
            fprintf('You must either correct your computations so that \n')
            fprintf('the difference d2c-d2c'' is of the order of 1E-16 at most \n')
            fprintf('or make d2c symmetric by doing the following trick: \n')
            fprintf('d2c = 0.5*(d2c+d2c''); \n')
            error('Cannot proceed with more testing until fixed.');
        end
    end
end

if PriLev == 1
    fprintf('No major problems detected. \n');
    fprintf('Do not call this function in production mode. \n');
end

exitFlag = 0;

% MODIFICATION LOG
%
% 050213 med   File written.
% 050215 med   Checks on x_L, x_U, x_0, A, b_L and b_U added.
% 050215 med   Updated linear constraint message and f nargout check
% 060814 med   FUNCS used for callbacks instead