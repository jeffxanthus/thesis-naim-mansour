%function [X, O, F, Fpen, F00, Cc, C, nCon, nSample, ExDText, initRS] =...
%          expDesign(Percent, nSample, AddMP, nTrial, CLHMethod, SCALE,...
%          RandState,PriLev,Prob,varargin)
%
% expDesign computes an initial experimental design for CGO solvers.
%
% INPUT PARAMETERS - Note! No error checking made, all input must be correct
%
% NEW CGO PARAMETERS:  nTrial, CLHMethod
%
% Percent     Type of strategy to get the initial sampled values:
%
%              Percent |  Exp. Design                   |  ExD
%             ---------|--------------------------------|-------
%                      |  CORNER STRATEGIES:            |
%                  900 |  All Corners                   |    1
%                  997 |  x_L + x_U + adjacent corners  |    2
%                  998 |  x_U + adjacent corners        |    3
%                  999 |  x_L + adjacent corners        |    4
%                      |                                |
%                      |  DETERMINISTIC STRATEGIES:     |
%                    0 |  User given initial points     |    5
%                99-94 |  Use glcDirect                 |    6
%                      |                                |
%                      |  LATIN BASED SAMPLING:         |
%                    1 |  Maximin LHD 1-norm            |    7
%                    2 |  Maximin LHD 2-norm            |    8
%                    3 |  Maximin LHD Inf-norm          |    9
%                    4 |  Minimal Audze-Eglais          |   10
%                    5 |  Minimax LHD (only 2 dim)      |   11
%                    6 |  Latin Hypercube               |   12
%                    7 |  Orthogonal Samling            |   13
%                      |                                |
%                      |  RANDOM STRATEGIES: (pp in %)  |
%                  1pp |  Circle surrounding            |   14
%                  2pp |  Ellipsoid surrounding         |   15
%                  3pp |  Rectangle surrounding         |   16
%
%             Negative values of Percent result in Constrained versions
%             of the Exp. Design methods 7-16. It means that all points
%             sampled are feasible with respect to all given constraints.
%
%             For ExD 5,6-12,14-16 user defined points are used
%
% nSample     Number of sample points to be used in initial experimental
%             design. nSample is used differently dependent on Percent:
%
%                     (n)Sample:
%               ExD |     < 0      |    = 0     |    > 0     |      []
%             ------|--------------|------------|------------|-------------
%                 1 |     2^d      |            |            |
%                 6 | abs(n) iters |            |            |
%              7-11 |     d+1      |    d+1     | max{d+1,n} | (d+1)(d+2)/2
%                12 |   LATIN(k)   |            |            |
%                13 |    abs(n)    |            |            |
%             14-16 |     d+1      |            |            |
%             ------|--------------|------------|------------|-------------
%
%             where LATIN = [21 21 33 41 51 65 65] and k = abs(nSample).
%
%             Otherwise nSample as input does not matter.
%
%             DESCRIPTION OF THE EXPERIMENTAL DESIGNS:
%
%             ExD 1, All Corners. Initial points is the corner points of the
%             box given by Prob.x_L and Prob.x_U. Generates 2^d points, which
%             results in too many points when the dimension is high.
%
%             ExD 2, Lower and Upper Corner point + adjacent points.
%             Initial points are 2*d + 2 corners: the lower left corner x_L
%             and its d adjacent corners x_L+(x_U(i)-x_L(i))*e_i, i=1,...,d
%             and the upper right corner x_U and its d adjacent corners
%             x_U - (x_U(i)-x_L(i))*e_i, i=1,...,d
%
%             ExD 3. Initial points are the upper right corner x_U and its
%             d adjacent corners x_U - (x_U(i)-x_L(i))*e_i, i=1,...,d
%
%             ExD 4. Initial points are the lower left corner x_L and its
%             d adjacent corners x_L + (x_U(i)-x_L(i))*e_i, i=1,...,d
%
%             ExD 5. User given initial points, given as a matrix in CGO.X.
%             Each column is one sampled point. If d = length(Prob.x_L),
%             then size(X,1) = d, size(X,2) >= d+1. CGO.F should be defined
%             as empty, or contain a vector of corresponding f(x) values.
%             Any CGO.F value set as NaN will be computed by solver routine.
%
%             ExD 6. Use determinstic global optimization methods to find the
%             initial design. Current methods available (all DIRECT methods).
%
%             Percent:   99 = glcDirect,    97 = glcSolve,    95 = glcFast
%                        98 = glbDirect,    96 = glbSolve,    94 = glbFast
%
%             ExD 7-11. Optimal Latin Hypercube Designs (LHD) with respect to
%             different norms. The following norms and designs are available:
%
%             Percent:    1 = Maximin 1-Norm,         4 = Audze-Eglais Norm
%                         2 = Maximin 2-Norm,         5 = Minimax 2-Norm
%                         3 = Maximin Inf-Norm,
%
%             All designs taken from:  http://www.spacefillingdesigns.nl/
%
%             Constrained versions will try bigger and bigger designs up to 
%             M = max{10*d,nTrial} different designs, stopping when it has
%             found nSample feasible points.
%
%             ExD 12. Latin hypercube space-filling design. For nSample < 0,
%             k = abs(nSample) should in principle be the problem dimension.
%             The number of points sampled is:
%                   k      :  2  3  4  5  6  >6
%                   Points :  21 33 41 51 65 65
%             The call made is: X = daceInit(abs(nSample),Prob.x_L,Prob.x_U);
%
%             Set nSample = []  to get (d+1)*(d+2)/2 sampled points:
%                   d      : 1  2  3  4  5  6  7  8  9 10
%                   Points : 3  6 10 15 21 28 36 45 55 66
%             This is a more efficient number of points to use.
%
%             If CGO.X is nonempty, these points are verified as in ExD 5,
%             and treated as already sampled points. Then nSample additional
%             points are sampled, restricted to be close to the given points.
%
%             Constrained version of Latin hypercube only keep points that
%             fulfill the linear and nonlinear constraints. The algorithm
%             will try up to M = max{10*d,nTrial} points, stopping when it
%             has found nSample feasible points (d+1 points if nSample < 0).
%
%             For pure integer programming (IP) problems with 1 or more 
%             linear equalities, a special function mipFeasible generates 
%             part of the design solving a set of LP problems. At least as 
%             many infeasible points as number of linear equalities are 
%             sampled, which is necessary.
%
%             ExD 13. Orthogonal Sampling, LH with subspace density demands.
%
%             ExD 14,15 and 16. Random strategies, the abs(Percent) value gives
%             the percentage size of an ellipsoid, circle or rectangle around
%             the so far sampled points that new points are not allowed in.
%             Range 1%-50%. Recommended values 10% - 20%.
%
%             If CGO.X is nonempty, these points are verified as in ExD 5,
%             and treated as already sampled points. Then nSample additional
%             points are sampled, restricted to be close to the given points.
%
% AddMP       If = 1, add the midpoint as extra point in the corner strategies.
%             Default AddMP=1 if any corner strategy.
%
% nTrial      For CLH, the method generates M = max{10*d,nTrial} trial points,
%             and evaluate them until nSample feasible points are found.
%
%             In the random designs, nTrial is the maximum number of trial
%             points randomly generated for each new point to sample.
%
% CLHMethod   Different search strategies for finding feasible LH points.
%             First of all, the least infeasible point is added. Then the
%             linear feasible points are considered. If more points are
%             needed still, the nonlinear infeasible points are added.
%
%             1 - Take the sampled infeasible points in order.
%             2 - Take a random sample of the infeasible points.
%             3 - Use points with lowest cErr.
%
% SCALE       0 - Original search space
%             1 - Transform search space to unit cube
%
% RandState   If >=0,  rand('state',RandState) is set to initialize the
%                      pseudo-random generator
%             If < 0   rand('state',sum(100*clock)) is set to give a new
%             or = []  set of random values each run.
%             If isnan(RandState), the random state is not initialized.
%
% PriLev      Print Level, if > 0 print 1 line information, if >1 more info
%
% Prob        Structure, where the following variables are used:
%   x_L       Lower bounds for each element in x.
%   x_U       Upper bounds for each element in x.
%   b_L       Lower bounds for the linear constraints
%   b_U       Upper bounds for the linear constraints
%   A         Linear constraint matrix
%   c_L       Lower bounds for the nonlinear constraints
%   c_U       Upper bounds for the nonlinear constraints
% --------------------------------------------
% optParam    Structure in Prob, Prob.optParam
% ---------------------------------------------
%             Defines optimization parameters. Fields used:
%  bTol       Linear constraint tolerance
%  cTol       Nonlinear constraint tolerance
%
% ------------------
% Fields in Prob.CGO
% ------------------
%             Defines user given points. Could also contain function values.
%             If ExD == 12,14-16 these points are included into the design.
%
% X           A matrix of initial x values. One column for every x value. 
%             If ExD == 5, size(X,2) >= dim(x)+1 needed.
% F           A vector of initial f(x) values. If any element is set
%             to NaN it will be computed.
% CX          Optionally a matrix of nonlinear constraint c(x) values.
%             If nonempty, then size(CX,2) == size(X,2). If any element
%             is set as NaN, the vector c(x) = CX(:,i) will be recomputed.
%
% varargin    Additional parameters that are sent to the costly f(x)
%
% ---------------------------------------
% MIP         Structure in Prob, Prob.MIP
% ---------------------------------------
%             Defines integer optimization parameters. Fields used:
%   IntVars:  Index vector with indices for integer variables
%
% OUTPUT PARAMETERS
%
% X           Matrix with sampled points (in unit space if SCALE == 1)
% O           Matrix with sampled points (in original space)
% F           Vector with function values (penalty added for costly Cc(x))
% Fpen        Vector with function values + additional penalty if infeasible
%             using the linear constraints and noncostly nonlinear c(x) 
% F00         Vector of pure function values, before penalties
% Cc          Matrix with costly constraint values, Cc(x)
% C           Matrix with noncostly constraint values, c(x)
% nCon        Number of noncostly constraint evaluations calling c(x)
% nSample     Corrected value of nSample
% ExDText     String with experimental design information
% initRS      State of random number generator immediately after initialization
%             If no initialization is done, its current value is saved instead.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Feb 6, 2006. Last modified Sept 19, 2009.

function [X, O, F, Fpen, F00, Cc, C, nCon, nSample, ExDText, initRS] =...
   expDesign(Percent, nSample, AddMP, nTrial, CLHMethod, SCALE,...
   RandState,PriLev,Prob,varargin)

if nargin < 9
   error('expDesign needs 9 parameters');
end

d         = Prob.N;
x_L       = Prob.x_L;
x_U       = Prob.x_U;
dLin      = Prob.mLin;
if Prob.simType == 1
   dCon   = 0;
else
   dCon   = Prob.mNonLin;
end

if isfield(Prob.MIP,'IntVars')
    IntVars   = Prob.MIP.IntVars;
else
    Prob.MIP.IntVars = [];
    IntVars = [];
end

bTol = Prob.optParam.bTol;      % Linear constraint feasibility tolerance
cTol = Prob.optParam.cTol;      % Constraint feasibility tolerance

x_D     = x_U - x_L;
% Scaling parameters
if SCALE
   x_LL = zeros(d,1);
   x_UU = ones(d,1);
   x_DD = x_UU;
else
   x_LL = x_L;
   x_UU = x_U;
   x_DD = x_D;
end

%NHQ    There exist new method for matlab versions 7.4 and later:
%       Mersenne Twister algorithm, rand('twister',RandState).
%       Stick to 'state' for now, maybe no need to change this??

% Set initial state for pseudo random generator
if isnan(RandState)
   % Do no initialization
elseif isempty(RandState)
   rand('state',sum(100*clock));
elseif (length(RandState)==1 && RandState >= 0) || length(RandState) == 35
   rand('state',RandState);
else
   rand('state',sum(100*clock));
end
% Vector with the current state of rand, used in the ExD. process.
% Even though RandState is Nan or [], we save the current value.
initRS = rand('state');

% Initialize output variables.
X = [];    F = [];    Fpen = [];
O = [];    C = [];    nCon = 0;
CX = [];  % Initial set of nonlinear constraint values

ExD  = [];
if ~isempty(ExD)
   %ExD = 0;
elseif isempty(Percent)
   % Set default strategy, Constrained Audze-Eglais LHD
   ExD = 10;
   Percent = -4;
elseif Percent == 900
   % Corner Strategy: Pick all corners of the box bounds
   ExD = 1;
elseif Percent == 997
   % Adjacent Corner Strategy: Pick x_L and x_U + all adjacent corner points
   ExD = 2;
elseif Percent == 998
   % Adjacent Corner Strategy: Pick x_U + all adjacent corner points
   ExD = 3;
elseif Percent == 999
   % Adjacent Corner Strategy: Pick x_L + all adjacent corner points
   ExD = 4;
elseif Percent == 0
   % Deterministic Designs: User given initial set of points
   ExD = 5;
elseif Percent > 90  &&  Percent < 100
   % Deterministic Designs: DIRECT
   ExD = 6;
elseif abs(Percent) == 1
   % Latin Based: Maximin Latin Hypercube Design, 1-norm
   ExD = 7;
elseif abs(Percent) == 2
   % Latin Based: Maximin Latin Hypercube Design, 2-norm
   ExD = 8;
elseif abs(Percent) == 3
   % Latin Based: Maximin Latin Hypercube Design, Inf-norm
   ExD = 9;
elseif abs(Percent) == 4
   % Latin Based: Minimax Latin Hypercube (only for 2 dim at the moment)
   ExD = 10;
elseif abs(Percent) == 5
   % Latin Based: Minimal Audze-Eglais Latin Hypercube
   ExD = 11;
elseif abs(Percent) == 6
   % Latin Based: Latin Hypercube and Constrained LH
   ExD = 12;
elseif abs(Percent) == 7
   % Latin Based: Orthogonal Sampling
   ExD = 13;
elseif abs(Percent) > 100  &&  abs(Percent) <= 150
   % Random Designs: Circle and Constrained Circle
   ExD = 14;
elseif abs(Percent) > 200  &&  abs(Percent) <= 250
   % Random Designs: Ellipsoid and Constrained Ellipsoid
   ExD = 15;
elseif abs(Percent) > 300  &&  abs(Percent) <= 350
   % Random Designs: Rectangle and Constrained Rectangle
   ExD = 16;
else
   fprintf('Invalid value of experimental design parameter Percent %d\n',Percent)
   error('expDesign: Illegal input!!!');
end


% Set the correct number of points to sample
%
%          (n)Sample:
%    ExD |     < 0      |    = 0     |    > 0     |      []
%  ------|--------------|------------|------------|-------------
%      1 |     2^d      |            |            |
%      6 | abs(n) iters |            |            |
%   7-11 |     d+1      |    d+1     | max{d+1,n} | (d+1)(d+2)/2
%     12 |   LATIN(k)   |            |            |
%     13 |    abs(n)    |            |            |
%  14-16 |     d+1      |            |            |
%  ------|--------------|------------|------------|-------------
%
%  where LATIN = [21 21 33 41 51 65 65] and k = abs(nSample).
%
%  Otherwise nSample as input does not matter.

nSample = nSampleSetup(ExD,d,nSample,AddMP);


% Go to the chosen Experimental Design
switch ExD
   case {1 2 3 4}  % All Corners + Double, Lower or Upper Adjacent Corners.
      if ExD == 1
         % Corner Strategy: Pick all corners of the box 
         X = corners(x_LL,x_UU);
         n = size(X,2);
         if ~isempty(nSample)
            % Use nSample many of the corners, select randomly
            [xxxx ix] = sort(rand(1,size(X,2)));
            X         = X(:,ix(1:nSample));
         else
            nSample = n;
         end
         if nSample == n
            ExDText = sprintf('ExD %d: All %d corners',ExD,n);
         else
            ExDText = sprintf('ExD %d: All Corners, Use %d of %d',ExD,nSample,n);
         end
      elseif ExD == 2
         % Adjacent Corner Strategy: Pick x_L and x_U + all adjacent corner points
         X = [adjacentL(x_LL,x_UU), adjacentR(x_LL,x_UU) ];
         if size(X,1) == 2
            % For d==2, only four corners
            X = X(:,1:4);
         elseif size(X,1) == 1
            % For d==1, only two corners
            X = X(:,1:2);
         end
         ExDText = sprintf('ExD %d: Double Adjacent Corners, n=%d',ExD,size(X,2));
      elseif ExD == 3
         % Adjacent Corner Strategy: Pick x_U + all adjacent corner points
         X = adjacentR(x_LL,x_UU);
         ExDText = sprintf('ExD %d: Upper Adjacent Corners',ExD);
      else
         % Adjacent Corner Strategy: Pick x_L + all adjacent corner points
         X = adjacentL(x_LL,x_UU);
         ExDText = sprintf('ExD %d: Lower Adjacent Corners',ExD);
      end

      % AddMP, this part identical for ExD 1,2,3 and 4.
      if AddMP
         OK = 1;
         % Avoid mid point if all variables are integer
         if length(IntVars) < d  ||  ( length(IntVars) == d  &&  sum(IntVars) < d )
            X = [X , 0.5*(x_LL +x_UU)];
            if ~isempty(IntVars)
               n = size(X,2);
               % Use experimental design values with equal probability integers
               X = daceInts(X,x_D,x_U,IntVars);
               if n > size(X,2)
                  OK = 0;
                  if PriLev > 1
                     fprintf('expDesign: Not possible to add midpoint.\n')
                  end
               end
            end
         end
      else
         OK = 0;
      end
      if OK
         ExDText = [ExDText,sprintf(' + midpoint.')];
      else
         ExDText = [ExDText,sprintf('.')];
      end
      
      % Scale back to original space
      if SCALE
         O = tomsol(9, x_L, X ,x_D);
      else
         O = X;
      end

      % Set dimension d and number of points n
      [d n]  = size(X);
      
      % Calculate function values. 
      [F00,Cc] = cgo_fc(O, Prob, varargin{:});
      F        = cgo_pf(F00, Cc, Prob);
      
      % HKH This could be parallellized
      if AddMP
         % midpoint might not have the lowest function value.
         if ~all( F(n) <= F(1:n-1) )  ||  OK == 0
            % generate "inner box corners"
            if ExD == 3
                xAdd = adjacentL(x_LL+0.25*x_DD,x_UU-0.25*x_DD);
            elseif ExD == 4
                xAdd = adjacentR(x_LL+0.25*x_DD,x_UU-0.25*x_DD);
            else
                xAdd = corners(x_LL+0.25*x_DD,x_UU-0.25*x_DD);
            end
            %NHQ must check here too if SCALE is set.
            if SCALE
               oAdd = tomsol(9, x_L, xAdd ,x_D);
            else
               oAdd = xAdd;
            end
            % Add a new point until interior point with lowest f
            nk = n;
            for k = 1:size(xAdd,2)
               nk      = nk+1;
               X(:,nk) = xAdd(:,k);
               O(:,nk) = oAdd(:,k);
               OK      = 1;
               if ~isempty(IntVars)
                  %X(:,nk) = sampleInts(X(:,nk),x_LL,x_UU,IntVars);
                  X = daceInts(X,x_D,x_U,IntVars);
                  % Check if "rounded xAdd" is already sampled
                  if nk > size(X,2)
                     nk = nk-1;
                     O  = O(:,1:nk);
                     OK = 0;
                  else
                     O  = X;
                  end
               end
               % Calculate function values and costly Cc. 
               if OK
                  [F00(nk),CcNew]    = cgo_fc(oAdd(:,k), Prob, varargin{:});
                  if isempty(CcNew)
                     F(nk)           = cgo_pf(F00(nk),[],Prob);
                  else
                     Cc(:,nk)        = CcNew;
                     F(nk)           = cgo_pf(F00(nk),CcNew,Prob);
                  end
                  % F(n+k) = nlp_f(oAdd(:,k), Prob, varargin{:});
               end
               if all( F(end) <= F(1:end-1) )
                  break
               end
            end
            ExDText = [ExDText,sprintf(' Additional %d points added,\ntrying to sample ',nk-n)];
            ExDText = [ExDText,sprintf('an interior f_min. In all %d points sampled.\n',nk)];
         end
      end
      
      if PriLev > 0
         fprintf('expDesign - ');
         fprintf('%s\n',ExDText);
      end

      

      %%%%%%%%%%%%%%%%%%%%%%% --- NEW CASE --- %%%%%%%%%%%%%%%%%%%%%%%
   case 5
      % User given initial set of points
      if isfield(Prob.CGO,'X')
         [OK,X,F,CX] = checkUserGivenPoints(Prob,SCALE);
         if ~OK
            error('Illegal set of user given points.')
            %return
         end
      else
         fprintf('The field Prob.CGO.X not set.\n')
         error('Illegal input to CGO solver.')
      end
      ExDText = sprintf('ExD %d: %d user given initial points.',ExD,size(X,2));
      if PriLev > 0
         fprintf('expDesign - ');
         fprintf('%s\n',ExDText);
      end

      %%%%%%%%%%%%%%%%%%%%%%% --- NEW CASE --- %%%%%%%%%%%%%%%%%%%%%%%
   case 6
      % Deterministic Designs: Use DIRECT methods as initial design.
      % nSample < 0 gives the number of function evaluations as the
      % result of doing abs(nSample) DIRECT iterations.
      if isfield(Prob.optParam,'MaxFunc')
         MaxFunc = Prob.optParam.MaxFunc;
      else
         MaxFunc = 200;
      end
      % Use a DIRECT method
      if nSample < 0
         % Set number of iterations
         Prob.optParam.MaxIter = abs(nSample);
         % Use Prob.optParam.MaxFunc that is total maximum
         minFunc = d+1;
      else
         % Set minimum number of function values
         Prob.optParam.MaxIter = 10000;
         Prob.optParam.MaxFunc = nSample;
         minFunc = nSample;
      end
      nMin = 0;
      nPnt = 0;
      while nPnt < minFunc && nMin < 100
         if nMin == 1
            % Warm start with one iteration more until enough func values
            Prob.WarmStart = 1;
            Prob.optParam.MaxIter = 1;
         end
         if nMin > 0
            Prob.optParam.MaxFunc = MaxFunc-nPnt;
         end
         switch Percent
            case 99
               GO = 'glcDirect';
               r=tomRunFast('glcDirect',Prob);
               warm=r.glcDirect.WarmStartInfo;
            case 98
               GO = 'glbDirect';
               r=tomRunFast('glbDirect',Prob);
               warm=r.glbDirect.WarmStartInfo;
               % Stupid name change until glbDirect is revised
               warm.C=warm.points;
               warm.F=warm.fPoints;
            case 97
               GO = 'glcSolve';
               tomRunFast('glcSolve',Prob);
               warm=load('glcSave.mat','C','F');
            case 96
               GO = 'glbSolve';
               tomRunFast('glbSolve',Prob);
               warm=load('glbSave.mat','C','F');
            case 95
               GO = 'glcFast';
               tomRunFast('glcFast',Prob);
               warm=load('glcFastSave.mat','C','F');
            case 94
               GO = 'glbFast';
               tomRunFast('glbFast',Prob);
               warm=load('glbFastSave.mat','C','F');
         end
         nMin = nMin + 1;
         % NHQ Avoid counting infeasaible points where
         % no function value has been calculated.
         %nPnt = length(warm.F);
         %nPnt = warm.nFuncTot;
         pIdx = find(warm.F < 1E5 );
         nPnt = length(pIdx);
      end
      % NHQ This part must be possible to simplify.
      % First the scaled [0,1] values warm.C are transformed to the full
      % space. Then if SCALE is set, they are transformed back to [0,1].
      % Also, if SCALE is set, the full space values are not saved and
      % hence calculated again in the end of this file.
      
      % NHQ Load only points where function value has been calculated
%       X = tomsol(9,x_L,full(warm.C),x_D);
%       F = warm.F;
      X = tomsol(9,x_L,full(warm.C(:,pIdx)),x_D);
      F = warm.F(pIdx);
      if SCALE
         D = x_D;
         D(D==0) = 1;
         for i=1:size(X,2)
            x = X(:,i);
            x = max(x_L,min(x_U,x));
            X(:,i) = (x-x_L)./D;
         end
      end
      
      % Any user-given points to be added?
      if isfield(Prob.CGO,'X')
         % If not OK, XX = [], FF = [] and CX = [] is returned.
         % If SCALE = 1, XX is scaled when returned.
         [tmp,XX,FF,CX] = checkUserGivenPoints(Prob,SCALE,0);
      else
         XX = [];
         FF = [];
         CX = [];
      end
      
      if ~isempty(XX)
         sX = size(X);
         X = [X XX];
         [X,uInd] = unique(X','rows');
         X = X';
%          O = [];    % New points added but no values for O.
%                     % Need to calculate O for all X.
      else
          uInd = 1:size(X,2);
      end
      if ~isempty(FF)
         F = [F ; FF];
         F = F(uInd);
      else
         F = [F ; NaN(size(XX,2),1)];
         F = F(uInd);
      end
      if ~isempty(CX)
         CX = [NaN(sX) CX];
         CX = CX(:,uInd);
      end
      
      ExDText = sprintf('ExD %d: %s Deterministic GO.',ExD,GO);
      if PriLev > 0
         fprintf('expDesign - ');
         fprintf('%s\n',ExDText);
      end


      %%%%%%%%%%%%%%%%%%%%%%% --- NEW CASE --- %%%%%%%%%%%%%%%%%%%%%%%
   case {7 8 9 10 11}
      % Optimal Latin Hypercube Designs
      %   Norm     Type               dim    nSample
      %  ----------------------------------------------
      %     1   Maximin L1-norm         2     2-inf
      %                                 3     2-16
      %     2   Maximin L2-norm       2-4     2-300
      %                              5-10     2-100
      %     3   Maximin Inf-norm        2     2-inf
      %                                 3     2-17,27-28,64-65
      %                               4-8     2-10
      %     4   Audze-Eglais (AE)    2-10     2-100
      %     5   Minimax Designs         2     2-27
      
      NORM = abs(Percent);
      normText = {'L1-Maximin LHD','L2-Maximin LHD','L-inf Maximin LHD', ...
                  'AE-Maximin LHD','Minimax LHD'};
      
      if Percent > 0
         Constrained = 0;
         % Every combination is not possible
         if d == 2  ||  ( d <= 10  &&  nSample <= 100 )
            if ( d < 2   ||   nSample < 2   ||   ~ismember(NORM,1:5) )                               || ...
               ( NORM == 1    &&  ( d >  3  ||  ( d == 3  &&  nSample > 16 ) ) )                     || ...
               ( NORM == 2    &&  ( d > 10  ||  ( ismember(d,2:4)   &&  nSample > 300 ) ...
                                            ||  ( ismember(d,5:10)  &&  nSample > 100 ) ) )          || ...
               ( NORM == 3    &&  ( d >  8  ||  ( d == 3  &&  ~ismember(nSample,[2:17,27:28,64:65]) ) ...
                                            ||  ( ismember(d,4:8)   &&  nSample > 10 ) ) )           || ...
               ( NORM == 4    &&  ( d > 10  ||  nSample > 100 ) )                                    || ...
               ( NORM == 5    && ~( d == 2  &&  nSample <= 27 ) )
               
               if PriLev > 1
                  disp('This combination is not supported right now. Switch to AE-norm.')
               end
               NORM = 4;
            end
            
            if NORM == 2
               X = getMaxiMinLHDnorm2(d,nSample);
            elseif NORM == 4
               X = getLHDnormAE(d,nSample);
            else
               X = getMaxiMinLHDnorm135(Percent,d,nSample);
            end
            
            N = size(X,2);
            % Horrible bug, making X in original space, using x_L and x_U
            %%X = repmat(x_L,1,N) + X.*repmat(x_D,1,N)/(N-1);
            %%%X = repmat(x_L,1,N) + ((X-1) + rand(d,N)).*repmat(xD,1,N)/N;
            X = repmat(x_LL,1,N) + X.*repmat(x_DD,1,N)/(N-1);
            %X = repmat(x_LL,1,N) + ((X-1) + rand(d,N)).*repmat(x_DD,1,N)/N;
         else
            if PriLev > 1
               fprintf('Doesn''t have a Maximin design that big. Switch to Standard LH.\n')
            end
            [X,O,F,Fpen,C,nCon] = expDesign(6, nSample,[], nTrial,CLHMethod, ...
                                                  SCALE,NaN,PriLev,Prob,varargin);
         end
         ExDText = sprintf(['ExD %d: ' normText{NORM} ', n=%d.'],ExD,size(X,2));
         
      else
         Constrained = 1;
         if  d == 2  ||  ( d <= 10  &&  nSample <= 100 )
            if isempty(nTrial)
               M = max(4000,10*max(d+1,nSample));
            else
               M = max([10*d,10*nSample,nTrial]);
            end
            
            % Find feasible LHD points
            [O,NORM,nCon] = constrainedLHD(NORM,nSample,Prob,nTrial,nCon,PriLev);
            OK = 1;
            
            % Same mistake here, making X in original space. Now fixed.
            if SCALE
               X = (O - repmat(x_L,1,size(O,2)))./repmat(x_D,1,size(O,2));
            else
               X = O;
            end
            
            % Have we found enough feasible points?
            if size(X,2) < nSample  &&  isempty(IntVars)
               % Include X in the LHD, and sample new points not too close to X
               X = daceInit(M,min(100,ceil(M/2)),x_LL,x_UU,X,0.05*norm(x_DD));
               OK = 0;
            end
         else
            if PriLev > 1
               fprintf('Doesn''t have a Maximin design that big. Switch to Standard LH.\n')
            end
            [X,O,F,Fpen,F00,Cc,C,nCon] = expDesign(-6, nSample,[], nTrial, ...
                                      CLHMethod,SCALE,NaN,PriLev,Prob,varargin);
            OK = 1;
         end
      end
      if ~isempty(IntVars)  &&  Constrained == 0
         % Use experimental design values with equal probability integers
         X  = daceInts(X,x_DD,x_UU,IntVars);
      end
      
      % Any user-given points to be added?
      if isfield(Prob.CGO,'X')
         % If not OK, XX = [], F = [] and CX = [] is returned.
         % If SCALE = 1, XX is scaled when returned.
         [tmp,XX,F,CX] = checkUserGivenPoints(Prob,SCALE,0);
      else
         XX = [];
         F  = [];
         CX = [];
      end
      
      if ~isempty(XX)
         sX = size(X);
         X = [X XX];
         [X,uInd] = unique(X','rows');
         X = X';
         O = [];    % New points added but no values for O.
                    % Need to calculate O for all X.
      end
      if ~isempty(F)
         F = [NaN(sX(2),1) ; F];
         F = F(uInd);
      end
      if ~isempty(CX)
         CX = [NaN(sX) CX];
         CX = CX(:,uInd);
      end
      
      if Constrained
         if OK == 0
            [X,O,nCon] = findFeasiblePoints(X,Prob,M,CLHMethod,PriLev,nSample,SCALE,nCon);
         end
         if nCon == 0
            ExDText = sprintf(...
                      ['ExD %d: Constrained ' normText{NORM} ', obtained %d.'],...
                      ExD,size(X,2));
         else
            ExDText = sprintf(...
                      ['ExD %d: Constrained ' normText{NORM} ', use %d get %d.'],...
                      ExD,nCon,size(X,2));
         end
      end
      
      if PriLev > 0
         fprintf('expDesign - ');
         fprintf('%s\n',ExDText);
      end


      %%%%%%%%%%%%%%%%%%%%%%% --- NEW CASE --- %%%%%%%%%%%%%%%%%%%%%%%
%    case {12 13}
   case 12
      % Latin Hypercube or Constrained LH.
      % If CLH   nSample < 0  gives  d+1 points
      % else     nSample < 0  k = abs(nSample) determines # of points:
      %                       k      :  2  3  4  5  6  >6
      %                       Points :  21 33 41 51 65 65
      
      XX   = [];
      dist = [];
      % Any user-given points?
      if isfield(Prob.CGO,'X')
         % If not OK, XX = [], F = [] and CX = [] is returned.
         [OK,XX] = checkUserGivenPoints(Prob,SCALE,0);
         if OK
            dist = 0.05*norm(x_DD);
         end
      end
      k = size(XX,2);
      
      if Percent <0 & length(IntVars) == d ...
         & (dCon == 0 & dLin > 0 & any(Prob.b_L==Prob.b_U))
         % Use mipFeasible for pure IP with linear constraints
         Constrained = 1;
         % Special when only linear constraints
         if isempty(nTrial)
            M = max(4000,10*max(d+1,nSample));
         else
            M = max([10*d,10*nSample,nTrial]);
         end
         Method = 1;
         X = mipFeasible(Prob, M, nSample, Method, PriLev);
         %[X,uInd] = unique(X','rows');
         %X = X';
         [L U P Q] = lu(sparse([X' ones(size(X,2),1)]));
         [p tmp_p] = find(P);
         % [L U p q] = lu(sparse([X' ones(size(X,2),1)]), 'vector');
         %X = X(:,p(1:nSample));
         %X = X(:,p(1:min(nSample,length(p))));
         Xs = X(:,p);
         
         if(rank(X) < d)
             fprintf('Warning: Rank of X (%i) is less than d (%i)', rank(X), d);
             X   = Xs(:,1:min(nSample,length(p)));   
             rnk = size(X,2);
         else
             nS = min(nSample, length(p));
             X = zeros(d, nS);
             rnk = 0;
             i   = 1;
             j   = 1;
             while(rnk < d) & j <= nS
                X(:,i) = Xs(:,j);
                nrnk = rank(X(:,1:i));
                if(nrnk > rnk)
                    i = i + 1;
                    j = j + 1;
                    rnk = nrnk;
                else
                    j = j + 1;
                end
             end
             while(rnk >= d) & j <= nS
                rnk = rnk + 1;
                X(:,rnk) = Xs(:,j);
                j = j + 1;
             end
             X = X(:,1:rnk);
         end
         
	 %xprinti(p,'p:');
	 %p1=rank([X' ones(size(X,2),1)])
	 %s1=svd([X' ones(size(X,2),1)]);
	 %q1=s1(end)/s1(1)
	 %for i=1:size(X,2)-1
         %    dXsny = min(tomsol(30,X(:,i),X(:,i+1:end)));
	 %    if dXsny == 0
	 %       fprintf('%d dist %f',i,dXsny);
	 %       fprintf(' ERROR!!!!!');
	 %       fprintf('\n');
         %    end
         %end
         if rnk < nSample
	    fprintf('expDesign:   Used %d points from mipFeasible\n',rnk);
            %X = [ X,daceInit(M,ceil(M/2),x_LL,x_UU)];
            %NHQ Include X in the LHD, and sample new points not too close to X
            X = daceInit(M,min(100,ceil(M/2)),x_LL,x_UU,[XX X],0.05*norm(x_DD));
	    fprintf('expDesign:   Added %d trial points with LH\n',size(X,2)-rnk);
         else
            X = [XX X];
         end

      elseif Percent < 0 &  ( dCon > 0 | dLin > 0 )
         Constrained = 1;
         %M = max(d+1,abs(Percent));
         if isempty(nTrial)
            M = max(4000,10*max(d+1,nSample));
         else
            M = max([10*d,10*nSample,nTrial]);
         end
         %X = daceInit(M, 100, x_LL, x_UU);
         X = daceInit(M,min(100,ceil(M/2)),x_LL,x_UU,XX,dist);
      else
         Constrained = 0;
         X = daceInit(nSample, [], x_LL, x_UU,XX,dist);
         if nSample < 0
            ExDText = sprintf('ExD %d: Latin Hypercube, d=%d, n=%d',...
                               ExD,d,size(X,2)-k);
         else
            ExDText = sprintf('ExD %d: Latin Hypercube, n=%d',...
                               ExD,size(X,2)-k);
         end
         if k > 0
            ExDText = [ExDText sprintf('+%d given.',k)];
         else
            ExDText = [ExDText '.'];
         end
      end
      if ~isempty(IntVars)
         % Use experimental design values with equal probability integers
         X  = daceInts(X,x_DD,x_UU,IntVars);
      end
      if Constrained
         % Constrained Latin Hypercube
         
         % DO NOT check user given points for feasibility!!!
	 % Must have infeasible points if equalities to get full rank
         if 0  % Yes, check all points
            [X,O,nCon] = findFeasiblePoints(X,Prob,M,CLHMethod,PriLev,nSample+k,SCALE,nCon);
            IDX = ismember(XX',X','rows');
            k = sum(IDX);
            X(:,end-sum(~IDX)+1:end) = [];
         else  % No, check only points found by LH.
            [X,O,nCon] = findFeasiblePoints(X(:,k+1:end),Prob,M,CLHMethod,PriLev,nSample,SCALE,nCon);
            [IDX,LOC] = ismember(XX',X','rows');
            X(:,LOC(IDX)) = [];
            X = [XX X];
            if k > 0
               if SCALE
                  O = [tomsol(9, Prob.x_L, XX, Prob.x_U-Prob.x_L) O];
               else
                  O = [XX O];
               end
            end
         end
         if nCon == 0
            ExDText = sprintf(...
                      'ExD %d: Constrained Latin Hypercube, obtained %d',...
	               ExD,size(X,2)-k);
         else
            ExDText = sprintf(...
                      'ExD %d: Constrained Latin Hypercube, use %d get %d',...
	               ExD,nCon,size(X,2)-k);
         end
         if k > 0
            ExDText = [ExDText sprintf('+%d given.',k)];
         else
            ExDText = [ExDText '.'];
         end
      end
      if PriLev > 0
         fprintf('expDesign - ');
         fprintf('%s\n',ExDText);
      end

      %%%%%%%%%%%%%%%%%%%%%%% --- NEW CASE --- %%%%%%%%%%%%%%%%%%%%%%%
   case 13
      % Orthogonal Sampling: Advanced LH
      % If COS   nSample < 0  gives  d+1 points
      % else        d :=
         %     N |  2    3    4    5    6    7    8     9    10
         %    ---|-----------------------------------------------
         %s := 2 |  4    8   16   32   64  128  256   512  1024
         %     3 |  9   27   81  243  729 2187 6561 19683 59049
      if Percent < 0 & (dCon > 0 | dLin > 0)
         Constrained = 1;
         %M = max(d+1,abs(Percent));
         if isempty(nTrial)
            M = 10*max(d+1,nSample);
         else
            M = max([10*d,10*nSample,nTrial]);
         end
         %X = daceInit(M, 100, x_LL, x_UU);
         %X = daceInit(M,ceil(M/2),x_LL,x_UU);
         X = orthSamp(M, [], x_LL, x_UU,2,PriLev);
      else
         Constrained = 0;
         if nSample < 0
            %          d :=
            %     N |  2    3    4    5    6    7    8     9    10
            %    ---|-----------------------------------------------
            %s := 2 |  4    8   16   32   64  128  256   512  1024
            %     3 |  9   27   81  243  729 2187 6561 19683 59049

            s = min( max(2,abs(nSample)) , 3);
            n = 1;
            X = orthSamp(n, s, x_LL, x_UU,2,PriLev);
            ExDText = sprintf('ExD %d: Orthogonal Sampling, d=%d, n=%d.',...
                               ExD,d,size(X,2));
         else
            X = orthSamp(nSample, [], x_LL, x_UU,2,PriLev);
            ExDText = sprintf('ExD %d: Orthogonal Sampling, n=%d.',...
                               ExD,size(X,2));
         end
      end
      if ~isempty(IntVars)
         % Use experimental design values with equal probability integers
         X  = daceInts(X,x_DD,x_UU,IntVars);
      end
      if Constrained
         % Constrained Orthogonal Sampling, use same method as for CLH
         [X,O,nCon] = findFeasiblePoints(X,Prob,M,CLHMethod,PriLev,nSample,SCALE,nCon);
         if nCon == 0
            ExDText = sprintf(...
                      'ExD %d: Constrained Orthogonal Sampling, obtained %d.',...
	               ExD,size(X,2));
         else
            ExDText = sprintf(...
                      'ExD %d: Constrained Orthogonal Sampling, use %d get %d.',...
	               ExD,nCon,size(X,2));
         end
      end
      if PriLev > 0
         fprintf('expDesign - ');
         fprintf('%s\n',ExDText);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%% --- NEW CASE --- %%%%%%%%%%%%%%%%%%%%%%%
   case {14 15 16}   %NHQ  Constrained AND Integer is now working.
      
      % Random Designs: Circle and Constrained Circle
      % Random Designs: Ellipsoid and Constrained Ellipsoid
      % Random Designs: Rectangle and Constrained Rectangle
      if Percent < 0  &  (dCon > 0 | dLin > 0)
         ExDText = sprintf('ExD %d: Constrained ',ExD);
         Constrained = 1;
         Percent = abs(Percent);
      else
         ExDText = sprintf('ExD %d: ',ExD);
         Constrained = 0;
      end

      % Determine which design to use. Find correct value for Percent.
      design  = ((-1)^Constrained)*floor(Percent/100);
      Percent = max(1 , min(50,Percent-100*abs(design)) );
      if abs(design) == 1
         ExDText = [ExDText, 'Random Spheres '];
      elseif abs(design) == 2
         ExDText = [ExDText, 'Random Ellipsoids '];
      else    %elseif design == 3
         ExDText = [ExDText, 'Random Rectangles '];
      end
      ExDText = [ExDText, sprintf('(%d%%)',Percent) ];
      if Constrained == 1
         if isempty(nTrial)
            nTrial  = 1000;
         end
         ExDText = [ExDText, sprintf(' Sample %d/each,', nTrial)];
      end

      % Any user-given points?
      if isfield(Prob.CGO,'X')
         % If not OK, XX = [], F = [] and CX = [] is returned.
         [OK,XX,F,CX] = checkUserGivenPoints(Prob,SCALE,0);
      else
         XX = [];
         F  = [];
         CX = [];
      end
      k = size(XX,2);
      
      [X,nCon] = randomDesigns(design,x_LL,x_UU,Percent,nSample,[],XX, ...
                                             nTrial,Prob,nCon,SCALE,PriLev);
      
      %NHQ possible error here. When using unique, X is sorted.
      % This is a big mistake if F is not empty, since its values
      % are not sorted like X. Tried to fix this now, must very
      % it though.
      [X,uInd] = unique(X','rows');
      X = X';
      sX2  = size(X,2);
      sDX2 = sX2 - k;
      if isempty(F)
         F = NaN(sX2,1);
      else
         F = [F ; NaN(sDX2,1)];
         F = F(uInd);
      end
      if ~isempty(CX)
         CX = [CX NaN(d,sDX2)];
      end
      
      if sDX2 < nSample
         ExDText = [ExDText, sprintf(' Found %d need %d.',sDX2,nSample)];
      else
         ExDText = [ExDText, sprintf(' Found %d',nSample)];
         if k > 0
            ExDText = [ExDText sprintf('+%d given.',k)];
         else
            ExDText = [ExDText '.'];
         end
      end
      if sX2 < d+1
         fprintf('expDesign: ');
         fprintf('Haven''t found enough initial points, interpolation will fail.\n')
         fprintf('Call to Constrained Latin Hypercube, in order not to crash.\n')
         fprintf('All feasible points found here will be included in sample.\n')
         
         % Save already found X and F values.
         Prob.CGO.X  = X;
         Prob.CGO.F  = F;
         Prob.CGO.CX = CX;
         
         % Call to LH or Constrained LH.
         LHmethod = ((-1)^Constrained)*6;
         [X,O,F,Fpen,F00,Cc,C,nCon2] = expDesign(LHmethod, nSample-sDX2, ...
                           [],nTrial,CLHMethod,SCALE,NaN,PriLev,Prob,varargin);
         
         nCon = nCon + nCon2;
         
         if size(X,2) > nSample + sDX2
            % Pick the last nSample points, avoiding taking away a
            % point where the costly F has already been calculated.
            X = X(:,end-nSample+1:end);
            F = F(end-nSample+1:end);
         end
      end
         
      if PriLev > 0
         fprintf('expDesign - ');
         fprintf('%s\n',ExDText);
      end

   otherwise
      error('Undefined Experimental Design Strategy.')
      %ExDText = ' ';
      %return
end

if isempty(O)
   % Scale back to original space?
   if SCALE
      O = tomsol(9, x_L, X ,x_D);
   else
      O = X;
   end
end

% Set dimension d and number of points n
[d n] = size(X);

if ~isempty(Prob.FUNCS.f)
   % Compute objective function value in each undefined point x in O
   if isempty(F)
      [F00,Cc]           = cgo_fc(O, Prob, varargin{:});
      F                  = cgo_pf(F00, Cc, Prob);
   else
      ix                 = find(isnan(F));
      [F00(ix,1),CcNew]  = cgo_fc(O(:,ix), Prob, varargin{:});
      if isempty(CcNew)
         Cc(:,ix)        = NaN;
         F(ix)           = cgo_pf(F00(ix),[],Prob);
      else
         Cc(:,ix)        = CcNew;
         F(ix)           = cgo_pf(F00(ix),CcNew,Prob);
      end
   end

   Fpen      = F;

   % Compute constraint values for all points
   if dLin > 0 
      L = Prob.A*O;
      for i = 1:n
         Fpen(i) = Fpen(i)+sum(max(0,...
                           max(Prob.b_L-bTol-L(:,i),L(:,i)-bTol-Prob.b_U)));
      end
   end
   if dCon > 0
      if isempty(CX)
         CX   = cgo_c(O, Prob, varargin{:});
         nCon = nCon + size(O,2);
      else
         ix   = find(any(isnan(CX)));
         if ~isempty(ix)
             CX(:,ix) = cgo_c(O(:,ix), Prob, varargin{:});
             nCon     = nCon + length(ix);
         end
      end
      for i = 1:n
          Fpen(i) = Fpen(i)+sum(max(0,...
                            max(Prob.c_L-cTol-CX(:,i),CX(:,i)-cTol-Prob.c_U)));
      end
      C = CX;
   end
else
   F00 = [];
   Cc  = [];
end

% ====================================================================
function nSample = nSampleSetup(ExD,d,nSample,AddMP)
% ====================================================================
%
%          (n)Sample:
%    ExD |     < 0      |    = 0     |    > 0     |      []
%  ------|--------------|------------|------------|-------------
%      1 |     2^d      |            |            |
%      6 | abs(n) iters |            |            |
%   7-11 |     d+1      |    d+1     | max{d+1,n} | (d+1)(d+2)/2
%     12 |   LATIN(k)   |            |            |
%     13 |    abs(n)    |            |            |
%  14-16 |     d+1      |            |            |
%  ------|--------------|------------|------------|-------------
%
%  where LATIN = [21 21 33 41 51 65 65] and k = abs(nSample).
%
%  Otherwise nSample as input does not matter.

switch ExD
   case {1}
      % Corner Strategy: Both x_L and x_U + adjacent corners
      if isempty(nSample)
         nSample = (d+1)*(d+2)/2;
      elseif nSample > 0
         nSample = max(d+1,nSample);
      elseif nSample == 0
         nSample = d+1;
      elseif nSample < 0
         nSample = 2^d;
      end
      nSample = min(2^d,nSample);

   case {2}
      % Corner Strategy: Both x_L and x_U + adjacent corners
      nSample = 2*d + 2 + AddMP;

   case {3 4}
      % Corner Strategy: x_L or x_U + adjacent corners
      nSample =   d + 1 + AddMP;

   otherwise
      % ExD = 6-11
      if isempty(nSample)
         nSample = (d+1)*(d+2)/2;
      elseif nSample > 0
         nSample = max(d+1,nSample);
      elseif nSample == 0
         nSample = d+1;
      elseif nSample < 0
         if ismember(ExD,[7:11,14:16])
            nSample = d+1;
         end
         %     elseif ExD == 6   glcDirect
         %         nSample = abs(nSample) iterations
         %     elseif ExD == 12  LH and CLH
         %         nSample = LATIN(k) or d+1
         %     elseif ExD == 13  OS and COS
         %         nSample = ORTH(k) or d+1
      end
end


% ====================================================================
function [OK,X,F,CX] = checkUserGivenPoints(Prob,SCALE,nCheck)
% ====================================================================
% This function checks that the matrix X of given points
% is valid in all ways necessary.
X  = Prob.CGO.X;
F  = [];
CX = [];

if nargin < 3
   nCheck = 1;
end

% Initialize flag to FALSE
OK = 0;

x_L = Prob.x_L;
x_U = Prob.x_U;
x_D = x_U - x_L;

d = Prob.N;
n = size(X,2);
dCon = Prob.mNonLin;
flag = [size(X,1) ~= d  n < d+1];
if nCheck == 0
   flag(2) = 0;
end
if any(flag)
   fprintf('Size of user given X matrix is wrong, size %d %d\n',size(X));
   fprintf('1st dimension must be   %d\n',d);
   fprintf('2nd dimension must be >= %d\n',d+1);
   if nCheck
      error('Illegal input to CGO solver')
   else
      X = [];
      return
   end
end
if isfield(Prob.CGO,'F')
   F = Prob.CGO.F(:);
end
if ~isempty(F) && length(F) ~= n
   fprintf('Length of user given F vector is wrong, size %d\n',length(X));
   fprintf('If nonempty must have length %d\n',n);
   if nCheck
      error('Illegal input to CGO solver')
   else
      X = [];
      F = [];
      return
   end
end
[X,ix] = rmduplic(X,x_L,x_U);
if size(X,2) < d+1  && nCheck
   disp('User-given points')
   disp(X)
   fprintf('Not enough unique points in X given, only %d\n',size(X,2));
   fprintf('Must be  > %d\n',d+1);
   error('Illegal input to CGO solver')
   %fprintf('Illegal input to CGO solver\n')
   %return
end
if ~isempty(F)
   F = F(ix);
end
if SCALE
   D = x_D;
   D(D==0) = 1;
   for i=1:size(X,2)
      x = X(:,i);
      x = max(x_L,min(x_U,x));
      X(:,i) = (x-x_L)./D;
   end
end
if isfield(Prob.CGO,'CX')
   CX = Prob.CGO.CX;
   if ~isempty(CX)
      if size(CX,1) ~= dCon && size(CX,2) ~= n
         disp('size of CX')
         size(CX)
         if nCheck
            error('expDesign: Illegal size of Prob.CGO.CX')
         else
            X  = [];
            F  = [];
            CX = [];
            return
         end
      end
   end
end
% If all tests passed, set flag to TRUE
OK = 1;


% ====================================================================
function X = corners(x_L,x_U)
% ====================================================================

d  = length(x_L);

ix = x_L~=x_U;
iV = find(ix);
iF = find(~ix);
m  = length(iV);

n = 2^m;

X = zeros(d,n);
for i = 1:length(iF)
   j = iF(i);
   X(j,:)=x_L(j);
end

for j=1:m
   var = iV(j);
   for i=1:2^j
      l=n/2^j;
      if mod(i,2)==1
         X(var,l*(i-1)+1:l*i)=x_L(var);
      else
         X(var,l*(i-1)+1:l*i)=x_U(var);
      end
   end
end

% ====================================================================
function X = adjacentL(x_L,x_U)
% ====================================================================

iV = find(x_L~=x_U);
m  = length(iV);
n  = length(iV)+1;
X  = x_L*ones(1,n);

for i=1:m
   X(i,i+1) = x_U(iV(i));
end

% ====================================================================
function X = adjacentR(x_L,x_U)
% ====================================================================

iV = find(x_L~=x_U);
m  = length(iV);
n  = length(iV)+1;
X  = x_U*ones(1,n);

for i=1:m
   X(i,i+1) = x_L(iV(i));
end

% ====================================================================
function X = sampleInts(X,x_L,x_U,IntVars)
% ====================================================================
% Not used now, instead daceInts
n = size(X,2);
for i = 1:length(IntVars)
   j = IntVars(i);
   X(j,:) = x_L(j) + floor(rand(1,n)*(x_U(j)-x_L(j)+1));
end
% Remove duplicates
ix = ones(n,1);
for i = 2:n
   if sum(all(X(:,1:i-1)==X(:,i)*ones(1,i-1))) > 0
      ix(i) = 0; % Mark duplicate point
   end
end
ix = find(ix);
if length(ix) < n
   % Only use unique points
   X = X(:,ix);
end

% ====================================================================
function [X,ix] = rmduplic(X,x_L,x_U)
% ====================================================================
% Remove duplicates
[d,n]  = size(X);
ix     = ones(n,1);
X(:,1) = max(x_L,min(x_U,X(:,1)));

for i = 2:n
   X(:,i) = max(x_L,min(x_U,X(:,i)));
   if any(sum(abs(X(:,1:i-1)-X(:,i)*ones(1,i-1)) < 1E-10)==d)
      ix(i) = 0; % Mark duplicate point
   end
end
ix = find(ix);
if length(ix) < n
   % Only use unique points
   X = X(:,ix);
end


% ====================================================================
function [X,O,nCon] = findFeasiblePoints(X,Prob,M,CLHMethod,PriLev,nSample,SCALE,nCon)
% ====================================================================
% Constrained Latin Hypercube. This function verifies and saves the
% feasible points. If not enough feasible points found, one of three
% possible strategies are used. Default is CLHMethod = 2.

% Default value of CLHMethod = 2
if isempty(CLHMethod)
   CLHMethod = 2;
elseif CLHMethod <= 0
   CLHMethod = 2;
end

dLin = Prob.mLin;
dCon = Prob.mNonLin;
bTol = Prob.optParam.bTol;      % Linear constraint feasibility tolerance
cTol = Prob.optParam.cTol;      % Constraint feasibility tolerance

%Scale back to original space before checking constraints
if SCALE
   O = tomsol(9, Prob.x_L, X, Prob.x_U-Prob.x_L);
else
   O = X;
end

nOK   = 0;
[d,n] = size(X);
ixOK  = zeros(n,1);

% Number of feasible PNTS to find
if nSample < 0
   PNTS = d+1;
else
   PNTS = nSample;
end

if PriLev > 0
   fprintf('Try <= %d to get %d feasible points, ',M,PNTS);
   fprintf('using CLHMethod %d.\n',CLHMethod);
end

if CLHMethod == 3
   cError = inf(n,1);
end
cErrBest = Inf;
ixBest   = Inf;

% loop and save feasible points
for i = 1:n
   if nOK >= PNTS
      break;
   end
   cErr = 0;
   if dLin > 0
      L = Prob.A*O(:,i);
      cErr = cErr+sum(max(0,max(Prob.b_L-bTol-L,L-bTol-Prob.b_U)));
   end
   if dCon > 0 && cErr == 0
      %C = nlp_c(O(:,i), Prob, varargin{:});
      C = nlp_c(O(:,i), Prob);
      nCon = nCon + 1;
      cErr = cErr+sum(max(0,max(Prob.c_L-cTol-C,C-cTol-Prob.c_U)));
   elseif dCon > 0
      ixOK(i) = -1;
   end
   if cErr == 0
      nOK = nOK + 1;
      ixOK(i) = 1;
   elseif cErr < cErrBest
      cErrBest = cErr;
      ixBest   = i;
   end
   if CLHMethod == 3
      cError(i) = cErr;
   end
end
% First add best infeasible point
if nOK < PNTS & ~isinf(ixBest)
   if ixOK(ixBest) == 0
      ixOK(ixBest) = 2;
      nOK          = nOK+1;
   end
end
if CLHMethod == 1
   % Take the first sampled infeasible points
   if nOK < PNTS
      % First use linear feasible points
      I       = find(ixOK == 0, PNTS-nOK);
      ixOK(I) = 3;
      nOK     = nOK+length(I);
   end
   if nOK < PNTS
      % Use linearly infeasible points
      I       = find(ixOK == -1, PNTS-nOK);
      ixOK(I) = 4;
   end
elseif CLHMethod == 2
   % Random sample of infeasible points
   if nOK < PNTS
      % First use linear feasible points
      I       = find(ixOK == 0);
      nI      = length(I);
      if nI <= PNTS-nOK
         ixOK(I) = 3;
         nOK     = nOK+nI;
      else % Random choice
         ix      = randperm(nI);
         I       = I(ix(1:PNTS-nOK));
         ixOK(I) = 3;
         nOK     = nOK+length(I);
      end
   end
   if nOK < PNTS
      % Use linearly infeasible points
      I       = find(ixOK == -1);
      nI      = length(I);
      ix      = randperm(nI);
      I       = I(ix(1:min(nI,PNTS-nOK)));
      ixOK(I) = 4;
   end
elseif CLHMethod == 3
   % Use points with lowest cErr
   if nOK < PNTS
      % First use linear feasible points
      I       = find(ixOK == 0);
      nI      = length(I);
      if nI <= PNTS-nOK
         % Use all linear feasible points
         ixOK(I) = 3;
         nOK     = nOK+nI;
      else % Choose points with lowest cErr
         [cE,ix] = sort(cError(I));
         I       = I(ix(1:PNTS-nOK));
         ixOK(I) = 3;
         nOK     = nOK+length(I);
      end
   end
   if nOK < PNTS
      % Use linearly infeasible points, choose points with lowest cErr
      I       = find(ixOK == -1);
      %nI      = length(I);
      [cE,ix]  = sort(cError(I));
      I       = I(ix(1:PNTS-nOK));
      ixOK(I) = 4;
   end
end

ix = find(ixOK > 0);
X  = X(:,ix);
O  = O(:,ix);

if PriLev > 0
   iV = find(ixOK == 1);
   fprintf('Found %d out of %d ',length(iV),PNTS);
   fprintf('feasible points at trials ');
   
   if isempty(iV)
      xprinti(iV,[],[],25)
   elseif iV(end) == PNTS
      fprintf('1 - %d.\n',PNTS);
   else
      xprinti(iV,[],[],25)
   end
   fprintf('\n');
   if length(iV) < PNTS
      fprintf('Deviation for best non feasible point %f\n',cErrBest);
   end
   if PriLev > 1
      fprintf('ixOK: ');
      fprintf('%d ',ixOK(ixOK > 0));
      fprintf('\n');
      if CLHMethod == 3
         xprint(cError(ixOK > 0),'cErr: ');
      end
   end
end

% MODIFICATION LOG
%
% 060206  hkh  First version, from code in arbfmip
% 060310  hkh  Call randEllips, new name for randomtest
% 060629  hkh  Add DIRECT methods as Percent 101 to 106
% 071004  hkh  Now 3 different adjacent corner strategies, -997,-998,-999
% 071006  hkh  Added SCALE and RandState as direct input parameters
% 071006  hkh  No random state init if isnan(RandState)
% 071008  hkh  Revised output for PriLev=1
% 080327  nhq  Major changes, restructuring the code.
%              New meaning of the parameters Percent and nSample.
%              New parameters nTrial and CLHMethod added.
%              New output parameter initRS.
%              New random design methods added.
% 080410  hkh  Revision of bTol, cTol. O output from findFeasiblePoints
% 080415  hkh  Wrong test for duplicates in sampleInts and
%              removed any point with 1 equal component
% 080415  hkh  Use error and not return in cases of error
% 080417  hkh  Add output string ExDText, used for result output
% 080424  hkh  Expanding mid point strategy
% 080531  hkh  max(4000,...) for CLH. Avoid crash for pure IP
% 080608  hkh  Revision of printing, ExDText for PrintResult
% 080613  hkh  Use neg Percent if both Prob.FUNCS.c and Prob.A empty
% 080614  hkh  New method: use LP/MIP solver to get feasible if linear problem
% 080614  hkh  Put priority on 1,2,...,d variables to obtain good points
% 080615  hkh  findFeasiblePoints: Must check if ixBest == inf
% 080625  hkh  Simplify output if all feasible points are in order 1,...,PNTS
% 080625  hkh  Use xprinti to print out order of feasible points found
% 080626  hkh  Use LU to find best points to use for pure IP problems
% 080629  hkh  Use cgo_fc to get vector of F values, possibly using >1 CPUS
% 080629  hkh  Handle Prob.simType = 0,1,2, costly/noncostly constraints
% 080630  hkh  Vectorize func/cons evals. Use cgo_c for noncostly constraints
% 080630  hkh  Add output Cc for costly constraints
% 080701  hkh  cgo_pf computes penalty function F, original F, F00, returned
% 080710  nhq  Added new designs, based on Maximin LHDs.
% 080710  nhq  Fixed Random Designs to handle Integers AND Constraints
% 080711  nhq  Now user given points with Random designs and LH-designs.
% 080717  hkh  Bug fixes, OK must be set if calling expDesign
% 080719  hkh  Must check for empty CcNew before update of Cc matrix
% 080804  hkh  Bug in new Percent 7-11, computing X in original user space
% 080916  nhq  Another bug in Percent 7-11, comp. X in original user space
% 080923  nhq  Found and Fixed bug in Percent 900 (an issue for MIP problems)
% 081007  nhq  Added the possibility to inlcude user given X ontop of
%              Maximin LHD designs. Uniqueness of the resulting X is fixed.
% 081008  nhq  Added the possibility to inlcude user given X ontop of
%              Global Solver. Uniqueness of the resulting X is fixed.
% 081101  nhq  When using DIRECT solvers, only count points where the
%              function value has actually been calculated. Fixed!!
% 081105  hkh  Use PriLev > 1 for "Switch to" printouts. 
% 081105  hkh  NORM = 4; must be set outside print statement
% 090424  nhq  Now possible to get design points even if Prob.FUNCS.f = []
% 090817  hkh  If isempty Prob.FUNCS.f avoid both f and constraint computation
% 090817  hkh  isempty Prob.FUNCS.f not possible for corner strategies ExD 1-4
% 090817  hkh  mipFeasible now only for Percent < 0, dCon == 0, dLin >0, 
%              length(IntVars) == d and 1 or more linear equalities
% 090901  hkh  Revised use of mipFeasible, rank determination and LH
% 090917  hkh  Limit # of intervals, s=min(100,M/2), in call to daceInit
% 090918  hkh  Revised daceInts generates equal probability integer values
