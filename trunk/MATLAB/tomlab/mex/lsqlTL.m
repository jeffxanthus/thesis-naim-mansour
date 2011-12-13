% TOMLAB LSQL special case of QP/LP Solver
%
% function Result = lsqlTL(Prob)
%
% INPUT:
% Prob   Problem structure in TOMLAB format
%
% x_L, x_U  Bounds on variables.
% b_L, b_U  Bounds on linear constraints.
% A         Linear constraint matrix.
% QP.c      Linear coefficients in objective function.
% QP.F      Quadratic matrix of size nnObj x nnObj.
% PriLevOpt Print level in solver (0-4)
%           0 - No output
%           1 - Only final error message
%
% -----------------------------------------------
% Fields used in Prob.LSQL:
% -----------------------------------------------
% PrintFile  Name of LSQL Print file. Amount and type of printing 
%            determined by PriLevOpt.
%
% OUTPUT:
% Result   Structure with results (see ResultDef.m):
% f_k      Function value at optimum.
% x_k      Solution vector.
% x_0      Initial solution vector.
% g_k      Exact gradient computed at optimum.
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
% ExitFlag Exit status from lsql.m (similar to TOMLAB).
% Inform   LSQL information parameter.
%           0 = Optimal solution with unique minimizer found
%           1 = Too many iterations
%           2 = Accuracy insufficient to attain convergence
%           3 = Internal inconsistency, division by zero
%           5 = An input parameter was invalid
%           else, constraint # not consistent with other active.
%           The problem has no feasible solution
%
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations. Set to Iter.
% GradEv   Number of gradient evaluations. Set to Iter.
% ConstrEv Number of constraint evaluations. Set to 0.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations. NOT SET.
% Solver   Name of the solver (lsql).
% SolverAlgorithm  Description of the solver.
%
% -----------------------------------------------------------------------
%
% For a problem description, see lsql.m
%
% -------------------------------------------------------------------------

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written Nov 1, 2000.  Last modified Aug 19, 2009.

function Result = lsqlTL(Prob)

if nargin < 1, error('lsqlTL needs the Prob structure as input');end

global MAX_x MAX_c % Max number of variables/constraints/resids to print

Prob.solvType = 2; % QP solver
Prob = iniSolveMini(Prob);

Result=ResultDef(Prob);
Result.Solver='LSQL';
Result.SolverAlgorithm='LSQL QP/LP code';

PriLev=Prob.PriLevOpt;

BIG=1E20;
[bl, bu, n, m] = defblbu(Prob, inf, 1);

% Initial checks on the inputs

x_0     = Prob.x_0(:);
if length(x_0) < n, x_0=zeros(n,1); end

% Safe-guard starting point
x_0   = max(bl(1:n),min(x_0,bu(1:n)));
x_L = bl(1:n);
x_U = bu(1:n);

% Make lower and upper bounds on linear constraints correct

if m > 0
   nA = size(Prob.A,2);
   if nA~=n, error('Linear constraints A MUST have n columns!'); end 
   b_L = bl(n+1:n+m);
   b_U = bu(n+1:n+m);
else
   b_L = [];
   b_U = [];
end

% Set up the constraint matrix A 

if m > 0
   ix1 = find( b_L == b_U & ~isinf(b_L) );
   ix2 = find( b_L ~= b_U & ~isinf(b_U) );
   ix3 = find( b_L ~= b_U & ~isinf(b_L) );
   
   A   = [full(Prob.A(ix1,:)); full(-Prob.A(ix2,:)); full(Prob.A(ix3,:))];
   b   = [-b_L(ix1); b_U(ix2); -b_L(ix3)];
   mEQ = length(ix1);
else
   A   = [];
   b   = [];
   mEQ = 0;
   ix1 = [];
   ix2 = [];
   ix3 = [];
end


if mEQ > 1
   % Indices of unbounded constraints
   xub_ix = find( x_L <= -BIG & x_U >= BIG);
   
   % Create a matrix B where the value of an element is a sum of the number
   % of nonzero elements in the corresponding row and column of every nonzero
   % element in A. Positions corresponding to a zero-value in A will have
   % the value Inf in B.
   % The lower the value, the more suitable the element in the corresponding
   % A-matrix is to be included in the identity matrix for the L special
   % constraints. 
   % Only the rows corresponding to equality constraints and the columns of 
   % unbounded variables are considered.
   AEQ = A(1:mEQ,xub_ix);
   AEQ_nz = find(AEQ);
   AEQ_ix = zeros(size(AEQ));
   AEQ_ix(AEQ_nz) = 1;
   temp = AEQ_ix*ones(size(AEQ,2))+ones(size(AEQ,1))*AEQ_ix;
   B = inf*ones(size(AEQ));
   B(AEQ_nz) = temp(AEQ_nz);

   % Build an array xy with the row- and column-indices of possible candidates
   [sorted_elements, ix] = sort(B(:));
   ix = ix(1:length(ix)-length(find(isinf(sorted_elements))));

   % Step through candidates and gradually increase the set composing the
   % eye-matrix.
   if ~isempty(ix)
      % Add first row
      xy = [mod(ix(1)-1,mEQ)+1,ceil(ix(1)/mEQ)];
      eyeix(1,:) = xy;
      norm_vec = AEQ(xy(1,1),xy(1,2));
      for i = 2:length(ix)
         % Get the row and column indices
         xy = [mod(ix(i)-1,mEQ)+1,ceil(ix(i)/mEQ)];
         include_flag = 1;
         if any(xy(1,1) == eyeix(:,1)) || any(xy(1,2) == eyeix(:,2))
            include_flag = 0;
         end
         for j = 1:size(eyeix,1)
            % Check if new row is free (has only zeros in previously included columns)
            if any(~isinf(B(xy(1,1),eyeix(:,2))))
               include_flag = 0;
            end
         end
         if include_flag == 1
            eyeix = [eyeix;xy];
            norm_vec = [norm_vec; AEQ(xy(1,1),xy(1,2))];
         end
      end
   else
      eyeix = [];
   end

   % The number of special equality constraints
   L = size(eyeix,1);

   % Only continue manipulating A and b if L > 1
   if L > 1
      % restore column indices (counting the bounded variables)
      eyeix(:,2) = xub_ix(eyeix(:,2));

      % Sort the rows and column and store the sorting order.
      % Also used for the result later.
      Arows = [1:size(A,1)]';
      Acols = [1:n]';
      row_sortorder = eyeix(:,1);
      for i = 1:size(A,1)
         if all(i ~= eyeix(:,1))
            row_sortorder = [row_sortorder; i];
         end
      end
      col_sortorder = [];
      for i = 1:n
         if all(i ~= eyeix(:,2))
            col_sortorder = [col_sortorder; i];
         end
      end
      col_sortorder = [col_sortorder; eyeix(:,2)];

      % Sort A, b, x_L and x_U
      A = A(row_sortorder,:);
      A = A(:,col_sortorder);
      b = b(row_sortorder);
      x_L = x_L(col_sortorder);
      x_U = x_U(col_sortorder);
      
      % Divide the rows with special constraints, and the corresponding bounds to make identity
      for i = 1:L
         Ai(i) = A(i,n-L+i);
         A(i,:) = A(i,:)/-Ai(i);
         b(i) = b(i)/-Ai(i);
      end
   else
      row_sortorder = 1:length(b);
      col_sortorder = 1:n;
   end
else
   L = 0;
   row_sortorder = 1:length(b);
   col_sortorder = 1:n;
end

H     = full(Prob.QP.F);
nnObj = size(H,1); % The number of nonlinear variables

% Check if any linear part
c = full(Prob.QP.c(:));

if L > 0 && nnObj > 0
   % Sort H and c
   H = H(col_sortorder,:);
   H = H(:,col_sortorder);
   c = c(col_sortorder);
end

% Determine type of problem
if isempty(c) | all(c==0)
   if isempty(H) | nnz(H) == 0
      Result.f_0=0;
   else
      Result.f_0=0.5*x_0(1:nnObj)'*H*x_0(1:nnObj);
   end
else
   if isempty(H) | nnz(H) == 0
      Result.f_0=c(1:n)'*x_0(1:n);
   else
      Result.f_0=0.5*x_0(1:nnObj)'*H*x_0(1:nnObj) + c(1:n)'*x_0(1:n);
   end
end

x_L(isinf(x_L)) = -BIG;
x_U(isinf(x_U)) =  BIG;

% Define default solver options.
LSQL      = DefPar(Prob,'LSQL',[]);
PriLevOpt = DefPar(Prob,'PriLevOpt',0);
PrintFile = DefPar(LSQL,'PrintFile','lsql.txt');

[x_k, Inform, Iter, iState, Ax, v, Obj] = lsql ( ...
        H, A, b, L, mEQ, c, x_L, x_U, PriLevOpt, PrintFile);

% lsql does not fix parameters on bounds, could be slightly off
x_k = min(x_U,max(x_L,x_k));

% Reorder multipliers to variables first
v = [v(size(A,1)+1:size(A,1)+n); v(1:size(A,1))];

% Restore the order of the solution variables
if L > 1
   x_k = x_k(sort(col_sortorder));
   v(1:n) = v(sort(col_sortorder));
   v(n+1:end) = v(n+sort(row_sortorder));
end

% Get the Lagrange multipliers back to the correct place in TOMLAB
if ~isempty(v)
   cLamda = zeros(n+m,1);
   cLamda(1:n) = v(1:n);
   cLamda(n+ix1) = v(n+1:n+length(ix1));
   j = n+length(ix1)+length(ix2);
   cLamda(n+ix2) = v(n+length(ix1)+1:j);
   for i = 1:length(ix3)
       if v(j+i) ~= 0
          cLamda(n+ix3(i)) = v(j+i);
       end
   end
end

if m > 0
   Ax = Prob.A*x_k;
end

switch Inform
   case 0
     ExitFlag=0;  % Convergence
   case 1
     ExitFlag=1;  % Too many iterations
   case {2,3}
     ExitFlag=3;  % Rank problem
   case 5
     ExitFlag=10; % Input errors
%  case -1
%    ExitFlag=2;  % Unbounded
   otherwise
     ExitFlag=4;  % Infeasible
end
if Inform == 2 & m > 0
   % Check if feasible
   if any(b_L - Prob.optParam.bTol > Ax) | any(b_U + Prob.optParam.bTol < Ax)
      ExitFlag = 4;
   end
end

Result.f_k = Obj;
Result.x_k = x_k;
Result.x_0 = x_0;
Result.v_k = cLamda;

if ~isempty(c)
   if ~isempty(H)
      Result.g_k=H*x_k+c;
   else
      Result.g_k=c;
   end
elseif isempty(H)
   Result.g_k=[];
else
   Result.g_k=H*x_k;
end

Result.FuncEv    = 0;
Result.GradEv    = 0;
Result.ConstrEv  = 0;
Result.Iter      = Iter;
Result.MinorIter = 0;
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

% Compute Result.xState, Result.bState and Result.QP.B only
Result = StateDef(Result, x_k, Ax, [], Prob.optParam.xTol, ...
         Prob.optParam.bTol, [], bl, bu, 1); 

switch Inform
   case 0
      Text = 'Optimal solution with unique minimizer found';
   case 1
      Text = 'Too many iterations';
   case 2
      Text = 'Accuracy insufficient to attain convergence';
   case 3
      Text = 'Internal inconsistency, division by zero';
   case 5
      Text = 'An input parameter was invalid';
   otherwise
      % If a constraint is conflictive, we need to map it to the original order
      % using row_sortorder
      if Inform-100 <= length(row_sortorder)
         Text = ...
         str2mat(['Constraint ' num2str(row_sortorder(Inform-100)) ...
                  ' not consistent with other active.'] ...
                 ,'The problem has no feasible solution');
      else
         Text = 'Unknown error';
      end
end

Result.ExitText = Text;

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nLSQL solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('LSQL: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
       fprintf('%s',Text(i,:))
       fprintf('\n')
   end
   fprintf('\n');
   
   fprintf('Objective function at solution x %36.18f\n\n',Obj);
   fprintf('Iterations      %7d. ',Iter);
   fprintf('\n');
   
   if PriLev > 1
      if isempty(MAX_x)
         MAX_x=n;
      end
      fprintf('Optimal x = \n');
      xprinte(x_k(1:min(n,MAX_x)),'x:  ');
   end
   if PriLev > 2
      fprintf('State vector iState for x and constraints = \n');
      xprinti(iState(1:min(length(iState),MAX_x)),'iState: ');
   end

   if PriLev > 3
      if isempty(MAX_c)
         MAX_c=20;
      end
      fprintf('Dual variables (Lagrangian multipliers) v_k (cLamda) = \n');
      xprinte(cLamda(1:min(length(cLamda),MAX_c)),'cLamda:');
   end
end
Result=endSolveMini(Prob,Result);

% MODIFICATION LOG:
%
% 090605 bjo  Written based on qlTL.m
% 090727 bjo  Added PrintFile
% 090819 bjo  Clean-up