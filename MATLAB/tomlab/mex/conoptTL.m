% TOMLAB CONOPT nonlinear constrained optimization solver
%
% function Result = conoptTL(Prob)
%
% CONOPT solves constrained nonlinear problems of the type:
%
%    min f(x)
%
%    subject to
%
%    x_L <= x <= x_U
%
%    g(x) <R> b
%
% where x,x_L,x_U are n-vectors and g(x) is an m-vector of linear and/or
% nonlinear functions.
%
% Furthermore, <R> designates exactly one of the relations:
%
%   <= (less than or equal to)
%   >= (greater than or equal to)
%   == (equal to)
%
% ---------------------------------------------------------------------
%
% The standard TOMLAB constrained nonlinear (con) problem definition is:
%
%    min f(x)
%
%    subject to
%
%    x_L <=  x   <= x_U, n  variable bounds
%    b_L <= A*x  <= b_U, m1 linear constraints
%    c_L <= c(x) <= c_U, m2 nonlinear constraints
%
% where x,x_L,x_U are n-vectors,
% b_L,b_U are m1-vectors, A is a m1xn matrix,
% c_L,c_L are m2-vectors and c(x) is an m2-vector of nonlinear functions.
%
% conoptTL and the interface functions cpt_c and cpt_dc automatically
% handle the translation from the TOMLAB formulation to the format needed
% by CONOPT, and back.
%
% This is done by the separation of constraints with both upper and lower bounds.
% For maximum performance, the user is advised to formulate his/her
% problems with single-bounded constraints if possible. Each double-bounded
% constraint is expanded to TWO distinct constraints internally.
%
% CONOPT is using 2nd order derivative information for the objective
% and the nonlinear constraints if any of the following is true:
% 1. Prob.NumDiff < 0 | Prob.ConsDiff < 0
% 2. Prob.FUNCS.H or Prob.FUNCS.d2c is nonempty, i.e. functions are given
%    to compute an analytic Hessian or 2nd order Lagrangian
%
% A special case is if the user has set Prob.CONOPT.LS2PTJ to 1.
% if 1 CONOPT is using internal 2nd order. In that case Tomlab is not
% providing or estimating any 2nd order information
%
% In the 2nd order case above CONOPT will estimate the pattern of the 2nd order
% Lagrangian, Prob.d2LPattern. It is more efficient if the user can provide
% this pattern
%
% -----------------------------------------------------------------------
%
% INPUT:
%
% Prob          Problem structure in TOMLAB format. Fields used are:
%
%  x_L, x_U     Bounds on variables.
%  b_L, b_U     Bounds on linear constraints.
%  A            Linear constraint matrix.
%  c_L, c_U     Bounds on nonlinear constraints.
%               For equality constraints (or fixed variables), set
%               e.g. b_L(k) == b_U(k).
%
%  ConsPattern  0-1 pattern telling the nonzero structure of the nonlinear
%               constraint Jacobian. One row for each nonlinear constraint.
%
%  d2LPattern   0-1 pattern telling the nonzero structure of the Hessian
%               of the Lagrangian function L(x,lam)=f(x)-sum(lam(i)*d2c_i(x))
%               If not given, and the # of variables is <300, a dense
%               Hessian is assumed. If more than 300 var's, second
%               derivatives are switched OFF.
%
%  PriLevOpt    Print level in solver.
%
%  optParam     Structure with optimization parameters. Many of these
%               parameters are communicated to the solver through
%               Prob.CONOPT.options, but ONLY if the corresponding field
%               in that structure is NOT already defined by the user.
%
%               Fields used:
%
%    MaxIter    Maximum number of iterations.
%               (Prob.CONOPT.options.LFITER)
%
%
%  CONOPT       Structure with special fields for CONOPT optimization
%               parameters. The following fields are used:
%
%    PrintFile  Name of file to print progress information and results to.
%               Default: '' (no file is opened).
%
%               To get progress information to the screen, set Prob.PriLevOpt = 1
%               (or higher).
%
%    StatFile   Name of file to print status information to.
%               Default: '' (no file is opened).
%
%    OptFile    Name of file with CONOPT options. Options set via this file
%               are overridden by options specified in Prob.CONOPT.options.
%               Default: Empty (no file used).
%
%    options    Structure with fields with names corresponding to
%               CONOPT options. The field names are case insensitive.
%
%               For example, to set the maximum number of iterations
%               to 10 000, do
%
%                  Prob.CONOPT.options.LFITER = 10000;
%
%               For a complete list of valid option names and their
%               meanings, see the users guide.
%
%    DebugFV    Set to nonzero to enable derivative debugging.
%               Default: 0 (off)
%
%               A more thorough but very expensive debugging mode can be
%               enabled by setting
%
%                  DebugFV = 1;
%
%               and also specifying
%
%                  Prob.CONOPT.options.LMDEBG = 1;
%
%    eqTol      A linear/nonlinear constraint is considered an equality
%               if the difference between the low and high limit is less
%               than this value. Default 1e-8.
%
% ----------------------------------------------------------------------------------
%
% OUTPUT:
% Result        Structure with optimization results.
%
%   f_k         Function value at optimum.
%   g_k         Gradient of the function.
%
%   x_k         Solution vector.
%   x_0         Initial solution vector.
%
%   c_k         Nonlinear constraint residuals.
%   cJac        Nonlinear constraint gradients.
%
%
%   xState      State of variables. Free == 0; On lower == 1; On upper == 2;
%               Fixed == 3;
%
%   bState      State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%               Equality == 3;
%
%   cState      State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2;
%               Equality == 3;
%
%   v_k         Lagrange multipliers (for bounds + dual solution vector).
%
%   ExitFlag    Exit status from CONOPT MEX.
%
%   Inform      CONOPT information parameter. See StatusText.
%
%   Iter        Number of iterations.
%   FuncEv      Number of function evaluations.
%   GradEv      Number of gradient evaluations.
%   ConstrEv    Number of constraint evaluations.
%
%   Solver           Name of the solver (CONOPT).
%   SolverAlgorithm  Description of the solver.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written May 20, 2003.   Last modified Sept 25 2011.

function Result = conoptTL(Prob)

%#function cpt_c cpt_dc cpt_d2L

if nargin < 1
    error('conoptTL needs the Prob structure');
end

Prob.solvType = 3; % NLP (CON) solver

CONOPT  = DefPar(Prob,'CONOPT',[]);
eqTol   = DefPar(CONOPT,'eqTol',1e-8);
options = DefPar(CONOPT,'options',[]);

LS2PTJ = 0;
if isempty(Prob.FUNCS.H) | (~isempty(Prob.FUNCS.c) & isempty(Prob.FUNCS.d2c))
   LS2PTJ = 1;
end
LS2PTJ  = DefPar(options,'LS2PTJ',LS2PTJ);

if LS2PTJ >= 1 | ~(Prob.NumDiff < 0 | Prob.ConsDiff < 0 | ~isempty(Prob.FUNCS.H) | ~isempty(Prob.FUNCS.d2c))
   % Internal estimation of 2nd order information by numerical differences
   Prob = iniSolve(Prob,3,1,1);
else
   % Use 2nd order information
   if ~isempty(Prob.FUNCS.H) & strmatch(funch2str(Prob.FUNCS.H),{'lp_H','qp_H'},'exact')
      if isempty(Prob.FUNCS.c) & isempty(Prob.d2LPattern)
         if ~isempty(Prob.QP.F)
            Prob.d2LPattern = spones(Prob.QP.F);
         else
            Prob.d2LPattern = sparse([], [], [], Prob.N, Prob.N, 0);
         end
         Prob = iniSolve(Prob,3,0,0);
      else
         if strmatch(funch2str(Prob.FUNCS.H),'qp_H','exact')
            if isempty(Prob.HessPattern)
               Prob.HessPattern = spones(Prob.QP.F);
            end
         else
            if isempty(Prob.HessPattern)
               Prob.HessPattern = sparse([], [], [], Prob.N, Prob.N, 0);
            end
         end
         Prob = iniSolve(Prob,3,2,2);
      end
   else
      Prob = iniSolve(Prob,3,2,2);
   end
end

Result = ResultDef(Prob);

% Algorithm text is set after call to solver
Result.Solver = 'CONOPT';

PriLev = DefPar(Prob,'PriLevOpt',0);

% Check bounds - although we'll disassemble lo,up again later
BIG = 1E6;
[lo,up,n,m1,m2] = defblbu(Prob,BIG,1);
m=m1+m2;

% Starting point
x_0 = DefPar(Prob,'x_0',zeros(n,1));

% Safe guard x_0
x_0 = max( lo(1:n), min( up(1:n),x_0)); 

Result.f_0 = nlp_f(x_0,Prob);
Result.x_0 = x_0;

A           = DefPar(Prob,'A',[]);
ConsPattern = DefPar(Prob,'ConsPattern',[]);

if isempty(ConsPattern)
   % No use complaining if m2 is zero 
   if Prob.Warning==1 & m2*n>1E6 % & (LargeScale | m2*n>=1E6)
      fprintf('\nWarning: ConsPattern is empty for largescale problem\nTo disable this message, set Prob.Warning=0\n\n');
   end
   ConsPattern = sparse(ones(m2,n));
end

if ~issparse(ConsPattern),ConsPattern = sparse(ConsPattern); end
if ~isempty(A) & ~issparse(A), A = sparse(A); end;

% Separate variable bounds and constraints
xl = lo(1:n);
xu = up(1:n);

bl = lo(n+1:n+m1);
bu = up(n+1:n+m1);

cl = lo(n+m1+1:end);
cu = up(n+m1+1:end);

% Linear constraint bounds
b     = []; btype   = [];
bex   = [];

for i=1:m1
   if abs(bl(i)-bu(i)) < eqTol       % Equality ?
      btype = [btype;0];
      b     = [b;bl(i)];
   elseif bl(i) > -BIG & bu(i) < BIG % Upper AND lower constrained?
      % Keep the upper bound
      btype = [btype;2];
      b     = [b;bu(i)];
      bex      = [bex;bl(i)]; % Add another constraint with the lower bound
      A        = [A;A(i,:)];  % Duplicate the current row, at the end of A
   elseif bu(i) >= BIG   % Upper bound larger than BIG - bounded below
      btype = [btype;1];
      b     = [b;bl(i)];
   elseif bl(i) <= -BIG  % Lower bound less than -BIG - bounded above
      btype = [btype;2];
      b     = [b;bu(i)];
   end
end

% All extra constraints are type 1 (lower bound)
b     = [b(:);bex(:)];
btype = [btype(:);ones(length(bex),1)];

% Nonlinear constraints
c     = []; ctype   = [];
cex   = [];

cexidx = 1:m2;
cexidx = cexidx(:);

for i=1:m2
   if abs(cl(i)-cu(i)) < eqTol  % Equality
      ctype  = [ctype;0];
      c      = [c;cl(i)];
   elseif cl(i) > -BIG & cu(i) < BIG  % Upper AND lower bounds
      % Keep the upper bound...
      ctype  = [ctype;2];
      c      = [c;cu(i)];
      % ... and add another constraint with the lower bound
      cex    = [cex;cl(i)];
      cexidx = [cexidx;i];
   elseif cu(i) >= BIG % Upper bound larger than BIG - bounded below
      ctype = [ctype;1];
      c     = [c;cl(i)];
   elseif cl(i) <= -BIG % Lower bound less than -BIG - bounded above
      ctype = [ctype;2];
      c     = [c;cu(i)];
   end
end

c     = [c;cex];
ctype = [ctype;ones(length(cex),1)]; % All extra constraints are type 1 (lower bound)

% Duplicate rows in ConsPattern for the added constraints
ConsPattern = [ ConsPattern ; ConsPattern(cexidx(m2+1:end),:) ];

Prob.m1     = m1;
Prob.m2     = m2;
Prob.cexidx = cexidx;

% The number of constraints in the extended problem
m1ex = m1+length(bex);
% m2ex = m2+length(cex);

% Avoid crash if user has wrongly set a c function, but no bounds
if m2 == 0 & ~isempty(Prob.FUNCS.c)
   Prob.FUNCS.c=[];
   Prob.FUNCS.dc=[];
   Prob.FUNCS.d2c=[];
end

nza = nnz(A);
nzc = nnz(ConsPattern);

if LS2PTJ < 1
   % Hessian information for Lagrangian
   HL = tril(DefPar(Prob,'d2LPattern',[]));

   % If no d2LPattern, estimate a dense hessian only if n<=1000. 
   % Otherwise switch off 2:nd derivatives. 
   if isempty(HL)
      if n<=1000
         HL = sparse(tril(ones(n,n))); 
         if(Prob.Warning)
            wstr = sprintf('\nWarning: No Prob.d2LPattern given. Using a full %d*%d Hessian pattern.',n,n);
            wstr = sprintf('%s\nTo disable this message, set Prob.Warning=0\n\n',wstr);
            fprintf(wstr);
         end
      else
         HL = [];
         if(Prob.Warning)
            wstr = sprintf('\nWarning: No Prob.d2LPattern given and n>=1000. 2nd derivatives disabled\n');
            wstr = sprintf('%s\nTo disable this message, set Prob.Warning=0\n\n',wstr);
            fprintf(wstr);
         end
      end
   end
   if ~issparse(HL), HL = sparse(HL); end
   % If there are ANY empty rows in the Hessian pattern, we have to fill
   % them out with ones - do this on the diagonal. 
   if any(sum(HL)==0)
      HL = HL + speye(size(HL));
   end
   nHess = nnz(HL);
else
   nHess = 0;
   HL    = [];
end

% Determine the sparse problem structure
if nHess > 0
   % Row and column indices. 
   % When used as arguments to the MEX, 1 should be subtracted.
   [hessrow,hesscol]=find(HL);
   
   % Send linear index from multiple subscripts for nonzero pattern
   Prob.D2Lidx = sub2ind(size(HL),hessrow,hesscol);
   hessrow = hessrow-1;
   hesscol = hesscol-1;
else
   Prob.D2Lidx = [];
   hessrow     = [];
   hesscol     = [];
end

PrintFile  = DefPar(CONOPT,'PrintFile','');
StatFile   = DefPar(CONOPT,'StatFile','');
OptFile    = DefPar(CONOPT,'OptFile','');
% options    = DefPar(CONOPT,'options',[]);
DebugFV    = DefPar(CONOPT,'DebugFV',0);

% Set some general Tomlab options in CONOPT.options - but only 
% if they're not already present. 

optParam = DefPar(Prob,'optParam',[]);
if isempty(optParam), optParam = optParamDef('CONOPT',Prob.probType,n,m); end

MaxIter = DefPar(optParam,'MaxIter',[]);
options.LFITER = DefPar(options,'LFITER',MaxIter);

MaxCPU = DefPar(Prob,'MaxCPU',1e6); % CONOPT Default value, 1E6
options.RVTIME = DefPar(options,'RVTIME',MaxCPU);

try
[modsta,solsta,iter,f_k,x_k,xmar,xbas,xsta,yval,ymar,ybas,ysta,vrsn] = ...
    conopt(...
    n,m1+length(bex),m2+length(cex),...
    x_0,xl,xu,...
    A,b,btype,nza,...
    c,ctype,ConsPattern,nzc,...
    nHess,hessrow,hesscol,...
    DebugFV,...
    PriLev,PrintFile,StatFile,OptFile,options,Prob) ;
catch
   l=lasterror;
   if(strcmp(l.identifier,'MATLAB:invalidMEXFile'))
      tomlabsharederror;
   else
      rethrow(l);
   end
end

Result.SolverAlgorithm = ['Feasible Path GRG, CONOPT ' vrsn];
Result.CONOPT.modsta = modsta;
Result.CONOPT.solsta = solsta;
Result.CONOPT.xval   = x_k;
Result.CONOPT.xmar   = xmar;
Result.CONOPT.xbas   = xbas;
Result.CONOPT.xsta   = xsta;
Result.CONOPT.yval   = yval;
Result.CONOPT.ymar   = ymar;
Result.CONOPT.ybas   = ybas;
Result.CONOPT.ysta   = ysta;

% Check for error
switch(modsta)
    case -999,
        
        % Something went very wrong - there is no solution available. 
        iter          = 0;
        
        % solsta is COI_Error in this case
        switch(solsta)
            case 105,   ExitText = 'Could not allocate basic memory needed to start CONOPT.';
            case 106,   ExitText = 'Too many equations. (Message contains the limit).';
            case 107,   ExitText = 'Too many variables. (Message contains the limit).';
            case 108,   ExitText = 'Too many equations + variables. (Message contains the limit).';
            case 109,   ExitText = 'Too many nonzeros. (Message contains the limit).';
            case 110,   ExitText = 'Model exceeds size limits for demonstration license. (Message contains the limits).';
            case 111,   ExitText = 'The index for the objective function defined with COIDEF_ObjVar or COIDEF_ObjCon is outside the legal range.';
            case 112,   ExitText = 'More nonlinear nonzeros than total nonzeros.';
            case 113,   ExitText = 'Insufficient memory to start CONOPT.';
            case 200,   ExitText = 'An internal system error has been encountered.';
            case 400,   ExitText = 'The model did not call Status and no other error has been seen. This may happen if the function or derivative debuggers find an error.';
                
                % These should never occur
            case 1002,  ExitText = 'Callback routine ReadMatrix was not registered.';
            case 1003,  ExitText = 'Callback routine FDEval was not registered.';
            case 1004,  ExitText = 'Callback routine ErrMsg was not registered.';
            case 1005,  ExitText = 'Callback routine Message was not registered.';
            case 1009,  ExitText = 'Callback routine Status was not registered.';
            case 1010,  ExitText = 'Callback routine Solution was not registered.';
        end
        
        % modsta values indicating "Normal" completion 
        
    case 1, ExitText = 'Optimal solution found';
    case 2, ExitText = 'Locally optimal';
    case 3, ExitText = 'Unbounded';
    case 4, ExitText = 'Infeasible';
    case 5, ExitText = 'Locally infeasible';
    case 6, ExitText = 'Intermediate infeasible';
    case 7, ExitText = 'Intermediate non-optimal';
        
        % case 8, ExitText = 'Not Used';
        % case 9, ExitText = 'Not Used';
        % case 10, ExitText = 'Not Used';
        % case 11, ExitText = '(Free for future use)';
        % case 14, ExitText = 'No Used';
        
    case 12, ExitText = 'Unknown type of error';
    case 13, ExitText = 'Error no solution';
    case 15, ExitText = 'Solved Unique';
    case 16, ExitText = 'Solved';
    case 17, ExitText = 'Solved Singular';
       
    otherwise, ExitText = 'Unknown modsta value';
end

% For all the "normal" modsta values, give some info about the solsta value
% too
if modsta == -999
    StatusText = 'Fatal error'; %, see Result.ExitText';
 else
    switch(solsta)
        case 1,  StatusText='Normal completion';
        case 2,  StatusText='Iteration interrupt';
        case 3,  StatusText='Resource interrupt';
        case 4,  StatusText='Terminated by solver';
        case 5,  StatusText='Evaluation error limit';
        case 6,  StatusText='Unknown';
        case 7,  StatusText='Not Used';
        case 8,  StatusText='User Interrupt';
        case 9,  StatusText='Error: Setup failure';
        case 10, StatusText='Error: Solver failure';
        case 11, StatusText='Error: Internal solver error';
        case 15, StatusText='Quick Mode Termination';
           
        otherwise, StatusText='Unknown SOLSTA value';
    end
end
Result.ExitText = [StatusText ' : ' ExitText];

if(modsta==1 | modsta==2)
   Result.ExitFlag = 0;
else
   Result.ExitFlag = solsta;
end

Result.Inform   = modsta;
Result.Iter = iter; 
Result.x_k  = x_k;
Result.f_k  = f_k;
Result.g_k  = nlp_g(x_k,Prob);
if m2>0
   % Result.c_k = yval(m1ex+1:m1ex+m2); DOES NOT GIVE CORRECT c_k
   Result.c_k  = nlp_c(x_k,Prob); 
   Result.cJac = nlp_dc(x_k,Prob);
else
   Result.c_k = [];
   Result.cJac = [];
end

% Lagrange multipliers, but for the original problem only. 
Result.v_k = [ xmar(1:n) ; ymar(1:m1) ; ymar(m1ex+1:m1ex+m2) ];

% ORIGINAL linear constraints at optimum
if m1>0 
    Result.Ax = yval(1:m1);
else
    Result.Ax = [];
end


Result = StateDef(Result, x_k, Result.Ax, Result.c_k, optParam.xTol, ...
                  optParam.bTol, optParam.cTol, lo, up, 1);

Result=endSolve(Prob,Result);

% MODIFICATION LOG
%
% 030520 ango Wrote file
% 030604 ango Robust handling of outputs added
% 030605 ango Revision of comments, options handling
% 030818 ango Minor changes
% 030901 ango Changed filename, used to be TconoptTL
% 031016 ango Require d2LPattern for second derivatives to be used, and
%             dense Hessian only if n<=300
% 031124 ango No Print/Status files are opened by default
% 031201 hkh  Avoid crash if c is given, but no lower or upper bounds
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040728 med  Pragmas added for MATLAB Compiler
% 041122 ango Updates for CONOPT 3.14, minor changes to warnings and info texts
% 041129 hkh  Use 1st or 2nd derivatives dependent on many things, see help
% 041130 ango Minor fix: Prob.FUNCS.d2c, not d2L
% 041201 ango Disable warning about missing ConsPattern for small problems
% 041202 hkh  Revise calls to defblbu, avoid unnecessary vector shuffling
% 050117 med  Options changed to options
% 050203 ango MaxCPU sent as RVTIME
% 050707 med  try-catch statement added to find install problem
% 050902 ango Handle problems with zero rows in Hessian
% 050908 med  Removed nlp_dc call unless m2 > 0
% 050908 med  d2LPattern set automatically for LP/QP/LPCON/QPCON problems before iniSolve call
% 050912 med  Added switch for LS2PTJ
% 050922 med  Updated installation error
% 050929 med  Removed xval and more
% 050929 med  CONOPT.eqTol added as control parameter
% 060217 med  Updated mex error catch
% 060814 med  FUNCS used for callbacks instead
% 070905 med  Fixed lp_H and qp_H checks for function handle
% 090804 med  Removed version check, cleaned
% 090925 hkh  Incorrect c_k computed, call nlp_c, and do it before calling nlp_dc
