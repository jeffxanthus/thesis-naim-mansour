% function stepexp(Name, t, y, wType, pInit, pEnd, SepAlg)
%
% Stepwise solve exponential fitting problems
%
% Name    Name of the problem
% t       Time steps t
% y       Observations y 
% wType   Weight type, 1 weight with data, 0 no weighting
% pInit   First p (number of exponential terms)
% pEnd    Last p  (number of exponential terms)
% SepAlg  If true, use separable nonlinear least squares
%
% stepexp calls both clsSolve and NLSSOL (if license is available)

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Apr 3, 2002.    Last modified Nov 23, 2004.

function Result = stepexp(Name, t, y, wType, pInit, pEnd, SepAlg)

% SepAlg = 1 implies the use of a separable least squares algorithm
% SepAlg  = 0; % Use ordinary least squares, no separation

fprintf('Run Time Series %s for %d terms to %d terms\n',Name,pInit,pEnd);

% Scale the problem
%yMax = max(y);
%y=y/yMax;    % Scale function values. Avoid large alpha. IMPORTANT!!!
%y = y / 100;

% Optimal values found:
x_opt=[];
f_opt =[];

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

% Select type of exponential model (see User's Guide)
eType=1;

% The following parameters are used by the initial value algorithm
% Normally do not change these values.

x0Type=0;
sumType=0;
infCR=0;
dType=0;
geoType=0;
qType=0;
sigType=0;

% Loop in the number of terms, p, from p0 to pEnd

fBest = inf*ones(pEnd,1);

fAlpha  = zeros(pEnd,pEnd);
fLambda = zeros(pEnd,pEnd);

for p=pInit:pEnd

fprintf('=============================================================\n');
fprintf('Run number of terms p = %d\n',p);
fprintf('=============================================================\n');
% Number of unknowns
n = 2*p;

% Just assign correct lengths for initial x_0, and bounds x_L and x_U
x_L = zeros(n,1);
x_U = ones(n,1);
x_0 = zeros(n,1);

x_min =x_L;
x_max =x_U;

Prob = clsAssign('exp_r', 'exp_J', JacPattern, x_L, x_U, Name, x_0, ...
                 y, t, weightType, weightY, SepAlg, fLowBnd, ...
                 A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ... 
                 x_min, x_max, f_opt, x_opt);

% Now reset the problem type, so Tomlab knows this is a exponential problem
% global probType
probType=checkType('exp');
Prob.probType=probType;     

Prob=expProbSet(Prob, p, wType, eType, infCR, dType, geoType,...
                qType, sigType, lambda, alpha, beta, x0Type, sumType);

ask = 0; 

% If ask = 1, interactive setting of:	
% SepAlg (ordinary versus separable nonlinear least squares strategy)
% wType  (no weighting, or weighting with data)
% p      (number of exponential terms)
% All initial values for lambda, one by one. Suggestion given by initial
% value algorithm in expInit.

% Find initial values
Prob=expInit(Prob,ask);

x_0 = Prob.x_0;

% Add linear constraints to robustify the optimization, avoiding
% singularities because of exponential components being too close. 

Prob = ExpLinCon(Prob);


[TomV,os,TV] = tomlabVersion;

if ~TV(3) | p < 3
   Res  = tomRun('clsSolve',Prob,2);

   fBest(p)       = Res.f_k;
   fLambda(1:p,p) = Res.x_k(1:p);
   fAlpha(1:p,p)  = Res.x_k(p+1:p+p);
   Result(p)      = Res;
end


if TV(3)
   if p < 3
      disp(' ')
      disp('Run nlssol')
      disp(' ')
      disp('Press return to continue')
      pause
   end
   % Here we may set parameters for the SOL solvers into the structure
   % See help nlssolTL for more help on the parameters
   Prob.SOL.optPar(1) = 11;  % Increase print level
   Prob.SOL.PrintFile = 'nlssol.txt';   % Name of text log print file
   Prob.SOL.SummFile  = 'nlssols.txt';  % Name of text log summary file

   % Special parameters for least squares problems are:
   % JTJ Initial Hessian (often best to have as true)
   Prob.SOL.optPar(47) = 1;  % Default 1, other unit Hessian

   % RESET Frequency . When Gauss-Newton works fine, often for small
   % residual problems, then one may raise this value
   Prob.SOL.optPar(48) = 2;  % Default 2, Reset each 2nd step

   R  = tomRun('nlssol',Prob,2);
   Res = R;
   if Res.f_k < fBest(p)
      fBest(p)       = Res.f_k;
      fLambda(1:p,p) = Res.x_k(1:p);
      fAlpha(1:p,p)  = Res.x_k(p+1:p+p);
      Result(p)      = Res;
   end
   if p < pEnd
      disp(' ')
      disp('Press return to continue')
      pause
   end
else
   disp(' ')
   disp('Press return to continue')
   pause
end


% Parameter statistics, get Jacobian from Result structure

Result(p).LS = StatLS(Result(p).x_k, Result(p).r_k, Result(p).J_k);

end

fprintf('Summary of Results\n')
fprintf('------------------\n')

N = length(y);

for p=pInit:pEnd
    fprintf('p = %d\n',p)
    fprintf('Weights =     ')
    fprintf('%12.5f ',fAlpha(1:p,p))
    fprintf('\n')
    fprintf('Intensities = ')
    fprintf('%12.5f ',fLambda(1:p,p))
    fprintf('\n')
end
fprintf('\n')
fprintf('Selection of Best Model\n')
fprintf('-----------------------\n')

for p=pInit:pEnd
    fprintf('p = %d ',p)
    fprintf('SSQ = %20.10f ',fBest(p)*2)
    fprintf('AIC = %12.5f ',AkaikeIC(fBest(p)*2,N,2*p))
    fprintf('BIC = %12.5f ',bic(fBest(p)*2,N,2*p))
    fprintf('AIC-corr = %12.5f ',aicc(fBest(p)*2,N,2*p))
    fprintf('\n')
end

