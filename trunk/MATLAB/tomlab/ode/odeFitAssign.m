% odeFitAssign is a direct way of setting up a parameter estimation
% problem in the TOMLAB (TQ) format.
%
% The information is put into the TOMLAB input problem structure Prob.
%
% function Prob = odeFitAssign(odeF, odeJ, Y0, E, Eeq, x_L, x_U, Name, ...
%                              t, x_0, odePattern, EWht, resWeight,...
%                              tInit, tStop, InitStep, ...
%                              A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ...
%                              x_min, x_max, f_opt, x_opt, ...
%                              IntVars, VarWeight, fIP, xIP);
%
% INPUT (Call with at least eight parameters)
%  odeF       Name of function that computes the ODE-system y'=f(y,t)
%             Option:
%             If odeF is a structure, assumed to be the Tomlab Prob structure
%             No call is then done to clsAssign, to define the Prob structure
%             The user must set Prob.ODE.f as odeF before or after the call
%             Only fields in Prob.LS and Prob.ODE are changed
%  odeJ       Name of function that computes the Jacobian of the ODE
%             system defined in odeF
%  Y0         Vector of length nEq (number of ODE equations) giving the
%             initial conditions of every ODE equation
%             If unknown, and to be estimated, set NaN
%             Tomlab stores a vector Prob.ODE.Y0Idx with pointers to these
%             unknowns. The first length(Y0Idx) elements in x are reserved
%             for these unknowns. In the input Prob structure to odeF:
%             Prob.ODE.Y0  The initial values Y0, with current x filled in
%             Prob.ODE.X   The other unknowns in x(length(Y0Idx)+1:end)
%  E          Matrix of observations. Each row represents a time value and
%             each column represents a series of measurements.
%             E(i,j) == NaN if an observation at time i is missing in
%             data series j.
%  Eeq        Index vector for matrix E. Eeq(i) = j means that
%             data series i is a measurement of ODE equation j, 1 <= j <= nEq
%             nEq = length(Y0)
%  x_L        Lower bounds on parameters. Default values is -inf
%  x_U        Upper bounds on parameters. Default values is inf
%
%      Note:  The number n of unknown parameters x is taken as
%             max{length(x_L), length(x_U), length(x_0)}.
%
%  Name       The name of the ODE parameter estimation problem.
% ........................................................
% The following parameters are optional, but recommended:
% ........................................................
%
%  t          Time values for the observations. length(t) must be the same
%             as the number of rows in E. t must be sorted in increasing order
%             Default if []: t = [1,2, ..., size(E,1)]
%
%  x_0        Initial values of parameters. Try to scale the problem so that
%             every x is roughly in the order of 1.
%             If [], x_0 = zeros(Default value is a vector of 0s.
%
%  odePattern 0-1 sparse or dense matrix that describes the
%             dependencies on x in each of the nEq (= length(Y0)
%             differential equations in odeF.
%             If odePattern(i,j) is 1, the equation i is dependent on x(j)
%
%  EWht       Weight on a complete series, length(EWht) == size(E,2)
%
%  resWeight  Weight on a single residual. The length of resWeight
%             must be the same as the number of NaN values in E.
%
%  tInit      Initial time Value. Default value is 0.
%
%  tStop      The ODE solver is not allowed to integrate past this point.
%             Default is empty.
%
%  InitStep   initial step length. Default set by the ODE solver
%
% L I N E A R   C O N S T R A I N T S
% A           Matrix A in linear constraints b_L<=A*x<=b_U. Dense or sparse.
% b_L         Lower bound vector in linear constraints, b_L<=A*x<=b_U.
% b_U         Upper bound vector in linear constraints, b_L<=A*x<=b_U.
%
% N O N L I N E A R   C O N S T R A I N T S
% c           Name of function that computes the mN nonlinear constraints
% dc          Name of function that computes the constraint Jacobian mN x n
% c_L         Lower bound vector in nonlinear constraints, c_L<=c(x)<=c_U.
% c_U         Upper bound vector in nonlinear constraints, c_L<=c(x)<=c_U.
% ConsPattern mN x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the constraint Jacobian and ones indicate values that
%             might be non-zero. Used when estimating the Jacobian numerically.
%             Estimated before solve, if Prob.LargeScale==1, ConsPattern==[]
%
% A D D I T I O N A L   P A R A M E T E R S
% x_min   Lower bounds on each x-variable, used for plotting
% x_max   Upper bounds on each x-variable, used for plotting
% f_opt   Optimal function value(s), if known (Stationary points)
% x_opt   The x-values corresponding to the given f_opt, if known.
%         If only one f_opt, give x_opt as a 1 by n vector
%         If several f_opt values, give x_opt as a length(f_opt) by n matrix
%         If adding one extra column n+1 in x_opt, 0 is min, 1 saddle, 2 is max.
%         x_opt and f_opt is used in printouts and plots.
%
% A D D I T I O N A L MIP P A R A M E T E R S
%
% If solving a mixed-integer nonlinear least squares problems
%
% IntVars     The set of integer variables. Can be given in one of three ways:
%
%             1) a scalar N<=n, in which case variables x(1)-x(N) are
%                restricted to integer values.
%
%             2) a vector of indices, e.g. [1 2 5]
%
%             3) a 0-1 vector of length n=length(x) where nonzero elements
%                indicate integer variables
%
% VarWeight   Priorities for each variable in the variable selection phase
%             A lower value gives higher priority.
%
% fIP         An upper bound on the IP value wanted. Makes it possible to
%             cut branches and avoid node computations.
% xIP         The x-values giving the fIP value.

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function Prob = odeFitAssign(odeF, odeJ, Y0, E, Eeq, x_L, x_U, Name, ...
                             t, x_0, odePattern, EWht, resWeight,...
                             tInit, tStop, InitStep, ...
                             A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ...
                             x_min, x_max, f_opt, x_opt, ...
                             IntVars, VarWeight, fIP, xIP)

if nargin < 32
   xIP=[];
   if nargin < 31
      fIP=[];
      if nargin < 30 
         VarWeight=[];
         if nargin < 29
            IntVars=[];
if nargin < 28
   x_opt=[];
   if nargin < 27
      f_opt=[];
      if nargin < 26
         x_max=[];
         if nargin < 25
            x_min=[];
            if nargin < 24
               c_U=[];
               if nargin < 23
                  c_L=[];
                  if nargin < 22
                     ConsPattern=[];
                     if nargin < 21
                        dc=[];
                        if nargin < 20
                           c=[];
                           if nargin < 19
                              b_U=[];
                              if nargin < 18
                                 b_L=[];
                                 if nargin < 17
                                    A=[];
end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end


if nargin < 16
   InitStep = [];
   if nargin < 15
      tStop = [];
      if nargin < 14
         tInit = [];
         if nargin < 13
            resWeight = [];
            if nargin < 12
               EWht = []
               if nargin < 11
                  odePattern = [];
                  if nargin < 10
                     x_0= [];
                     if nargin < 9
                        t = [];
                        if nargin < 8
                           error('Number of parameters must be >= 8')
end,end,end,end,end,end,end,end,end,

nPar = max([length(x_L), length(x_U), length(x_0)]);

if nPar < 1
  error('There must at least be one parameter to estimate');
end

% CALLCLS: 
% if == 0 Skip call to clsAssign, just define fields in Prob.LS, Prob.ODE

if isstruct(odeF)
   Prob = odeF;
   odeF = '';
   CALLCLS = 0;
else
   CALLCLS = 1;
end

nEq = length(Y0);

if nEq <  1
   error('Length of Y0 (number of ODE equations) must be > 0');
end

if isempty(E)
   error('Data matrix E must not be empty');
else
   E            = [NaN *ones(1, size(E,2)); E];
   [Em, nSerie] = size(E);
   if isempty(Eeq)
      error('Input parameter Eeq must not be empty');
   elseif length(Eeq) ~= nSerie
      fprintf('Input parameter Eeq must have length size(E,2)=%d\n',nSerie);
      error('Illegal length of input parameter Eeq');
   else
      Eeq = Eeq(:);
   end
   D            = zeros(Em,nSerie);
   for i=1:nSerie
       D(:,i) = (Eeq(i)-1)*Em + [1:Em]';
   end
   D            = D(:);
   E            = E(:);
   EIdx         = find(~isnan(E));
   y            = E(EIdx);
   yIdx         = D(EIdx);
   nRes         = length(EIdx);
end


if isempty(t)
   tWant = [0:Em];
   % error('Input parameter t must not be empty');
else
   if any(t~=sort(t))
      error('t is NOT an increasing sequence')
   else
      tWant = [0; t(:)]; % Must safeguard if t row vector
   end
end
T = repmat(tWant,1,nSerie);
T = T(:);
T = T(EIdx);

if isempty(EWht)
   % OK
elseif length(EWht) ~= nSerie
   fprintf('Length of EWht %d, should be %d\n',length(EWht),nSerie);
   error('Illegal length of EWht');
else
   EWht = EWht(:);
end

if isempty(resWeight)
   % OK
elseif size(resWeight,2) > 1 & size(resWeight,1) > 1
   if size(resWeight,1) ~= Em-1 | size(resWeight,2) ~= nSerie
      fprintf('resWeight: rows %d cols %d\n', size(resWeight)); 
      fprintf('E:         rows %d cols %d\n', Em-1,nSerie); 
      error('resWeight as matrix must have same size as E');
   end
elseif length(resWeight) ~= Em-1
   fprintf('Length of resWeight %d, should be %d\n', length(resWeight), Em-1); 
   error('Illegal length of resWeight');
else
   resWeight = resWeight(:); 
end

if isempty(resWeight) & (isempty(EWht) | all(EWht==1) )
   weightY    = []; 
   weightType = 0;
else
   weightType = 2;
   if size(resWeight,2) == 1
      resWeight = repmat([ 1;resWeight],1,nSerie);
   elseif isempty(resWeight)
      resWeight = ones(Em,nSerie);
   else
      % Add extra first row
      resWeight = [ones(1,nSerie);resWeight];
   end
   if ~isempty(EWht)
      resWeight = resWeight*diag(EWht);
   end
   resWeight = resWeight(:);
   weightY   = resWeight(EIdx);
end

JacPattern = [];
% Only make JacPattern if different series
% if nSerie > 1 & ~all(Eeq)==Eeq(1))
% end
%if ~isempty(odePattern)
%   JacPattern = zeros(nRes, nPar);
%   count = 1;
%   for i=1:nSerie
%      JacPattern(count :count + elE(i) - 1, :)  = ...
%      repMat(odePattern(Eeq(i),:), elE(i), 1);
%      count = count + elE(i);
%   end
%else
%   JacPattern = [];
%end

if isempty(tInit)
   tInit = 0;
end

if tInit > t(1)
   error('tInit must be less or equal to t(1)');
end

t(1) = tInit;

if tStop < max(t)
   error('tStop must be larger than or equal to max(t)');
end

if CALLCLS > 0
   % This must be changed if we add DAE and PDE 
   JacPattern = [];
   % Send a correct time vector T
   Prob = clsAssign('odeFit_r', [], JacPattern, x_L, x_U, Name, x_0, ...
                    y, T, weightType, weightY, [],0, ...
                    A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ...
                    x_min, x_max, f_opt, x_opt, ...
                    IntVars, VarWeight, fIP, xIP);
else
   % Just set the important fields in Prob.LS
   Prob.LS.y          = y;
   Prob.LS.t          = T;
   Prob.LS.weightType = weightType;
   Prob.LS.weightY    = weightY;
end

% Standard ODE parameters
if CALLCLS > 0
   Prob.ODE.f          = deblank(odeF);
   Prob.ODE.J          = deblank(odeJ);
elseif ~isempty(odeJ)
   Prob.ODE.J          = deblank(odeJ);
end
Prob.ODE.Y0         = Y0(:);
Prob.ODE.tWant      = tWant;       % Time points wanted
Prob.ODE.tInit      = tInit;
Prob.ODE.tStop      = tStop;
Prob.ODE.InitStep   = InitStep;
Prob.ODE.Name       = deblank(Name);

% Parameter estimation ODE parameters
Prob.ODE.Y0Idx      = find(isnan(Y0(:)));
Prob.ODE.nSeries    = nSerie;
Prob.ODE.odePattern = odePattern;
Prob.ODE.Eeq        = Eeq(:);
Prob.ODE.odePattern = odePattern;
Prob.ODE.yIdx       = yIdx;
Prob.ODE.EIdx       = EIdx;

% Needed for special case of CALLCLS == 0 and for MODFIT
Prob.ODE.t         = t;           % Time points corresponding to data E
Prob.ODE.E         = E;           % Data E
Prob.ODE.EWht      = EWht;
Prob.ODE.resWeight = resWeight;

Prob.probType      = checkType('ode');

% MODIFICATION LOG:
%
% 050413  joho Written
% 050417  hkh  Complete revision
% 050417  hkh  Add Y0 instead of nEq as input, Y0Idx NaN in Y0, nEq implicit
% 050418  hkh  Efficient new yIdx, old yIdx now EIdx
% 050418  hkh  weightType = 2 if using input weights
% 050418  hkh  Added linear and nonlinear constraints, and MIP information
% 050418  hkh  1st argument Prob, CALLCLS ==0, no clsAssign call
% 050421  bkh  Added check of tInit and tStop.
% 050422  bkh  Changed odeH_s/hStart, odeT_s/tStart and odeT_e/tEnd to
%              InitStep, tInit and tStop respectively
% 050503  bkh  Set probType to ODE
% 050705  med  Help updated