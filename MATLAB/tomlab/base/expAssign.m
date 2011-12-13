% function Prob = expAssign(p, Name, t, y, wType, eType, SepAlg, x_0)
%
% Create exponential fitting problems for given number of terms p
%
% Call Result = tomRun(Solver,Prob,2); to solve the problem.
%
% INPUT:
% p       Number of exponential terms
% Name    Name of the problem
% t       Time steps t
% y       Observations y, (t,y) must have the same length
% wType   Weight type, 1 = weight with data y, 0 =  no weighting
%
% The following parameters are optional:
%
% eType   Type of exponential model expression. Default 1, values 1-5.
%         (see User's Guide)
% SepAlg  If true, use separable nonlinear least squares. Default 0
% x_0     Initial values. If empty, TOMLAB initial value algorithm is used
%         x_0 = [lambda;alpha;beta]; where
%         lambda   p-dimensional lambda vector (intensities)
%         alpha    p-dimensional alpha  vector (weights)
%         beta     p-dimensional beta   vector (2nd weights when eType==4);
%         Note: always give full x_0, even if the alpha and beta values will
%         not be used at all when SepAlg = 1.
%
% OUTPUT:
% Prob    Exponential fitting problem
%
%      Notes:
%
%      Recommended setting for NLSSOL (are used by default when calling
%      expSolve).
%
%      Special parameters for least squares problems are:
%      JTJ Initial Hessian (often best to have as true)
%      Prob.SOL.optPar(47) = 1;  % Default 1, other unit Hessian
%
%      RESET Frequency . When Gauss-Newton works fine, often for small
%      residual problems, then one may raise this value
%      Prob.SOL.optPar(48) = 2;  % Default 2, Reset each 2nd step

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Aug 11, 2002.    Last modified Aug 28, 2006.

function Prob = expAssign(p, Name, t, y, wType, eType, SepAlg, x_0)

if nargin < 8
   x_0 = [];
   if nargin < 7
      SepAlg = [];
      if nargin < 6
         eType = [];
         if nargin < 5
            wType = 0;
            if nargin < 4
               error('At least 4 inputs required');
            end
         end
      end
   end
end

if isempty(SepAlg), SepAlg  = 0; end
if isempty(eType),  eType    = 1; end

eType =max(1,min(5,eType));

% SepAlg = 1  % Implies the use of a separable least squares algorithm
% SepAlg = 0; % Use ordinary least squares, no separation

%yMax = max(y);
%y=y/yMax;    % Scale function values. Avoid large alpha. IMPORTANT!!!
%y = y / 100;

% Optimal values ([] if not known, only used for print out):
x_opt = [];
f_opt = [];

% The routines exp_r and exp_J computes the exponential fitting residual
% and Jacobian for the given type of exponential model (eType)

% No special pattern is the Jacobian, it is a full matrix
JacPattern = [];

% For exponential problem, use routine ExpFitW for weighting
weightType=3;
weightY='ExpFitW';  % Define function to compute weights

% Lower bound on optimal function value
fLowBnd = 0;

% If linear constraints are present, set these in A, b_L and b_U
A = []; b_L = []; b_U = [];

% If nonlinear constraints are present, set these the bounds in c_L and c_U
c_L = []; c_U = [];

% The nonlinear constraint routine is c, contraint Jacobian is dc and
% the derivative pattern is set in ConsPattern
c = []; dc = []; ConsPattern = [];

% Set exponential parameters
lambda=[]; alpha=[]; beta=[];

% The following parameters are used by the initial value algorithm
% Normally do not change these values.

x0Type  = 0;
sumType = 0;
infCR   = 0;
dType   = 0;
geoType = 0;
qType   = 0;
sigType = 0;

% Number of unknowns
n = 2*p;

% Just assign correct lengths for initial x_0, and bounds x_L and x_U
x_L   = zeros(n,1);
x_U   = ones(n,1);
x_min = x_L;
x_max = x_U;

Prob = clsAssign('exp_r', 'exp_J', JacPattern, x_L, x_U, Name, x_0, ...
   y, t, weightType, weightY, SepAlg, fLowBnd, ...
   A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ...
   x_min, x_max, f_opt, x_opt);

% Now reset the problem type, so Tomlab knows this is a exponential problem
% global probType
probType      = checkType('exp');
Prob.probType = probType;

Prob=expProbSet(Prob, p, wType, eType, infCR, dType, geoType,...
   qType, sigType, lambda, alpha, beta, x0Type, sumType);

% If ask = 1, interactive setting of:
% SepAlg (ordinary versus separable nonlinear least squares strategy)
% wType  (no weighting, or weighting with data)
% p      (number of exponential terms)
% All initial values for lambda, one by one. Suggestion given by initial
% value algorithm in expInit.

% Find initial values
Prob = expInit(Prob);

if ~isempty(x_0)
   Prob.x_0 = full(double(x_0(:)));
end
if SepAlg
   Prob.x_0 = Prob.x_0(1:p);
   if ~isempty(Prob.A) & size(Prob.A,2) ~= Prob.N
      Prob.A = Prob.A(:,1:Prob.N);
   end
end

% Add linear constraints to robustify the optimization, avoiding
% singularities because of exponential components being too close.

Prob = ExpLinCon(Prob);

% MODIFICATION LOG:
%
% 041123  hkh  Change call to tomRun
% 050802  hkh  Add eType as input
% 050802  hkh  Handle SepAlg correctly
% 060130  med  Code separated to expAssign and expSolve
% 060822  med  All vectors set to full and double