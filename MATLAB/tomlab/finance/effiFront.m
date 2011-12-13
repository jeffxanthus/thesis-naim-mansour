% effiFront.m:
% 
% Computes the Mean Variance Efficient Frontier, or portfolios on the
% constrained efficient frontier
%
% function [pRisk, pRet, pWts] = effiFront(ExpRet, Cov, nPnt, RetVal, PLOT, ...
%                                x_L, x_U, A, b_L, b_U, Prob)
%
% To find Minimal Risk, a simplified QP is solved:
%
%        min  x' F x, where F is n x n,  F=Cov, x asset weights
%         x
%        s/t   x_L <=   x         <= x_U, x is in R^n
%                0 <= ExpRet(:)'x <= inf (constraint is never active)
%                1 <= ones(1,n) x <= 1
%              b_L <= A x         <= b_U, A has m rows
%              c_L <= c(x)        <= c_U  (NOT YET)
%
%        Minimal return is then ExpRet(:)' x (=minRet)
%        The portfolio variance is sqrt(x' F x )
%
%
% To find maximal possible return (maxRet), an LP problem is solved:
%
%        min  -ExpRet(:)' x  (= -maxRet)
%         x
%        s/t   x_L <=   x         <= x_U, x is in R^n
%                0 <= ExpRet(:)'x <= inf (constraint is never active)
%                1 <= ones(1,n) x <= 1
%              b_L <= A x         <= b_U, A has m rows
%              c_L <= c(x)        <= c_U  (NOT YET)
%
% Minimization problem: F=Cov, v = RetVal(i), i=2,NPnt-1 or i=1:length(RetVal)
%
%        min  x' F x, where F is n x n
%         x
%        s/t   x_L <=   x         <= x_U, x is in R^n
%                v <= ExpRet(:)'x <= v
%                1 <= ones(1,n) x <= 1
%              b_L <= A x         <= b_U, A has m rows
%              c_L <= c(x)        <= c_U  (NOT YET)
%
%        I.e. lowest risk given the expected return v, v in [minRet,maxRet]
%
% The model assumption is that asset returns are jointly normal.
% The portfolio return is then then ExpRet(:)' x 
% The portfolio variance is then sqrt(x' F x ) (=sqrt(x' Cov x ) )
% ---------------------------------------------------------------------------
%
% INPUT PARAMETERS
%   ExpRet     Expected Returns, length(ExpRet) = size(Cov,1) = size(Cov,2)
%   Cov        Estimated covariance matrix
%   nPnt       Number of points to compute on the efficient frontier,
%              default 20 (number of efficient portfolios)
%              Note! nPnt is only used if RetVal is []
%              nPnt equally spaced expected return values between the minimum 
%              and maximum values are computed (interval [minRet,maxRet])
%              If nPnt <= 2, only the minRet and maxRet portfolio is computed
%   RetVal     Vector with target portfolio return values on the frontier
%   PLOT       If > 0, plot efficient frontier
%   x_L        Vector with lower bounds on the asset weights x, default 0s
%   x_U        Vector with upper bounds on the asset weights x, default 1s
%   A          Matrix with linear constraints on the assets 
%              E.g. to specify an asset group, set 1:s if the asset is
%              included in the group, 0 otherwise
%   b_L        Vector with lower bounds on every linear asset constraint
%              if empty, assumed -inf
%   b_U        Vector with upper bounds on every linear asset constraint
%              if empty, assumed  inf
%   Prob       Optional TOMLAB problem structure, fields used:
%              SolverLP   LP solver
%              SolverQP   QP solver
%              PriLevOpt  Print level in tomRun calls
%
% OUTPUT PARAMETERS
%   pRisk      Vector with the standard deviation of return (risk) for each 
%              portfolio on the efficient frontier
%   pRet       Expected return for each portfolio
%   pWts       Matrix with the optimal weights for each portfolio.
%              Each column j holds the weights corresponding to portfolio j 
%
% EXAMPLE:
%
% See effiFrontEx.m in tomlab\finance

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Dec 18, 2004.   Last modified Dec 19, 2004.

function [pRisk, pRet, pWts] = effiFront(ExpRet, Cov, nPnt, RetVal, PLOT,...
                               x_L, x_U, A, b_L, b_U, Prob)

if nargin < 11
   Prob = [];
   if nargin < 10
      b_U = [];
      if nargin < 9
         b_L = [];
         if nargin < 8
            A = [];
            if nargin < 7
               x_U = [];
               if nargin < 6
                  x_L = [];
                  if nargin < 5
                     PLOT = [];
end, end, end, end, end, end, end

n       = size(Cov,1);
if n ~= size(Cov,2)
   error('The covariance matrix Cov is not quadratic');
end
if n ~= length(ExpRet)
   error('Rows/columns of covariance matrix Cov not same as length of ExpRet');
end
mG      = size(A,1);
PriLev  = DefPar(Prob,'PriLevOpt',[]);

if isempty(PriLev), PriLev = 1; end
if isempty(PLOT), PLOT = 1; end
if isempty(x_L), x_L = zeros(n,1); end
if isempty(x_U), x_U = ones(n,1); end
if mG > 0
   if isempty(b_L), b_L = -inf*ones(mG,1); end
   if isempty(b_U), b_U =  inf*ones(mG,1); end
end

A   = [ExpRet(:)';ones(1,n);A];
b_L = [0  ;1;b_L(:)];
b_U = [inf;1;b_U(:)];
x_L = x_L(:);
x_U = x_U(:);
F   = Cov + Cov';           % Make symmetric, now F = 2*Cov

ProbLP = lpAssign(-ExpRet, A, b_L, b_U, x_L, x_U, [], 'Max Return LP');

if ~isempty(RetVal)
   nPnt = length(RetVal);
else
   nPnt = max(2,nPnt);
end
pRisk = zeros(nPnt,1);
pRet  = zeros(nPnt,1);
pWts  = zeros(n,nPnt);

%             
SolverLP   = DefPar(Prob,'SolverLP',[]);
if isempty(SolverLP)
   SolverLP = GetSolver('lp', n > 1000, 1);
end
Result       = tomRun(SolverLP, ProbLP, PriLev);

ExitFlag     = Result.ExitFlag;
if ExitFlag == 1
   error('Error solving maxRet LP problem, problem is unbounded')
elseif ExitFlag == 4
   error('Error solving maxRet LP problem, problem is infeasible')
end

x_k          = Result.x_k;
pWts(:,nPnt) = x_k;
pRisk(nPnt)  = sqrt(0.5*(x_k'* F * x_k));
MaxRet       = -Result.f_k;          % MaxRet    = x_k'*ExpRet;
pRet(nPnt)   = MaxRet;

SolverQP     = DefPar(Prob,'SolverQP',[]);
if isempty(SolverQP)
   SolverQP  = GetSolver('qp', n > 300, 1);
end
Prob   = qpAssign(F, zeros(n,1), A, b_L, b_U, x_L, x_U, [], 'Mean Variance');
Prob   = ProbCheck(Prob,SolverQP,2,2);

Result       = tomRun(SolverQP, Prob, PriLev);
ExitFlag     = Result.ExitFlag;
if ExitFlag == 1
   error('Error solving MinRet QP problem, problem is unbounded')
elseif ExitFlag == 4
   error('Error solving MinRet QP problem, problem is infeasible')
end

x_k          = Result.x_k;
MVRet        = x_k'*ExpRet;
MVRisk       = sqrt(Result.f_k);

% Some printing, could be made optional
if 1
   fprintf('\n');
   fprintf('effiFront:       LP solver %s\n',SolverLP);
   fprintf('effiFront:       QP solver %s\n',SolverQP);
   fprintf('Minimal possible return is %f\n',MVRet)
   fprintf('Found at risk value        %f\n',MVRisk)
   fprintf('Maximal possible return is %f\n',MaxRet)
   fprintf('Found at risk value        %f\n',pRisk(nPnt))
   fprintf('\n');
end

if isempty(RetVal)
   if MVRet >= MaxRet
      nPnt   = 1;
      pWts   = pWts(:,nPnt);
      pRisk  = pRisk(nPnt);
      pRet   = pRet(nPnt);
      return
   else
      pWts(:,1)      = x_k;
      pRet(1)        = MVRet;
      pRisk(1)       = MVRisk;
      pRet(2:nPnt-1) = [MVRet+(1:nPnt-2)*(MaxRet-MVRet)/(nPnt-1)]';
      ix             = [2:nPnt-1];
   end
else
   if any(MaxRet < RetVal)
      fprintf('Minimal possible return is %f\n',MVRet)
      fprintf('Found at risk value        %f\n',pRisk(1))
      fprintf('Maximal possible return is %f\n',MaxRet)
      error('Some values in RetVal is bigger than Maximal possible return') 
   end
   if any(MVRet > RetVal)
      fprintf('Minimal possible return is %f\n',MVRet)
      fprintf('Found at risk value        %f\n',pRisk(1))
      fprintf('Maximal possible return is %f\n',MaxRet)
      error('Some values in RetVal is lower than minimal possible return') 
   end
   pRet = RetVal(:);
   ix   = [1:nPnt];
end

for i=ix
    Prob.x_0    = x_k;       % Use previous solution as initial point
    % No need to pose an equality constraint, upper bound may be inf
    Prob.b_U(1) = pRet(i); 
    Prob.b_L(1) = pRet(i); 
    Result      = tomRun(SolverQP, Prob, PriLev);
    ExitFlag    = Result.ExitFlag;
    if ExitFlag == 1
       fprintf('Failure, solving for Expected Return %f\n',pRet(i));
       error('Error solving QP problem, problem is unbounded')
    elseif ExitFlag == 4
       fprintf('Failure, solving for Expected Return %f\n',pRet(i));
       error('Error solving QP problem, problem is infeasible')
    end
    x_k         = Result.x_k;
    pWts(:,i)   = x_k;
    pRisk(i)    = sqrt(Result.f_k);
end

if PLOT > 0
   EffFront = figure;

   set(EffFront,'NumberTitle','off');
   set(EffFront,'Name','Efficient Frontier');
   set(EffFront,'Resize','on');
	
   set(EffFront,'Tag','EffFront');
	
   plot(pRisk, pRet);
   title('Mean-Variance Efficient Frontier', 'Color', 'k');
   xlabel('Risk(Standard Deviation)');
   ylabel('Expected Return');
   grid on;
      
   % prevent any output
   clear EffFront;
end

% MODIFICATION LOG:
%
% 021218  hkh  Written
