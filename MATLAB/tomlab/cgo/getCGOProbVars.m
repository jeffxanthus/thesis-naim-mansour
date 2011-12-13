%function [SCALE,REPLACE,RandState,Percent,nSample,AddMP,nTrial,CLHMethod,... 
%          SMOOTH,globalSolver,localSolver,LOCAL,GOMaxFunc,GOMaxIter,...
%          maxFunc1,maxFunc2,maxFunc3,GOlocalSolver,DIRECT,IntVars,Reals,...
%          nMax,d,x_L,x_U,x_D] = getCGOProbVars(Prob, x_L, x_U)
%
% getCGOProbVars picks up the input in Prob.CGO:
%
% Default values are set for
%   SCALE, RandState, SMOOTH, globalSolver, localSolver, GOlocalSolver, DIRECT
% but not for
%   REPLACE, Percent, nSample, AddMP, nTrial, CLHMethod
%
% getCGOProbVars sets the values of:
%   IntVars, Reals, d using Prob structure information
%
% x_L, x_U are safe guarded to integer values if MIP problem, x_D = x_U-x_L;
%
%
% If a pure IP problem, nMax is computed, otherwise nMax = Inf
%
% LOCAL - depends on SMOOTH and if IntVars is [] or not.
%
% The Prob.GO variables: GOMaxFunc,GOMaxIter,maxFunc1,maxFunc2,maxFunc3
% are extracted from the Prob structure, and set if pure IP problem.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written April 10, 2008.    Last modified May 31, 2008.

function [SCALE,REPLACE,RandState,Percent,nSample,AddMP,nTrial,CLHMethod,... 
          SMOOTH,globalSolver,localSolver,LOCAL,GOMaxFunc,GOMaxIter,...
          maxFunc1,maxFunc2,maxFunc3,GOlocalSolver,DIRECT,IntVars,Reals,...
          nMax,d,x_L,x_U,x_D] = getCGOProbVars(Prob, x_L, x_U)

x_D   = x_U - x_L;
if sum(x_D) == 0
   error('All variables are fixed')
end
d     = length(x_L);  % Problem dimension
nMax  = Inf;

if isempty(Prob.CGO)
   SCALE        = []; REPLACE      = []; RandState   = []; 
   Percent      = []; nSample      = []; AddMP       = []; 
   SMOOTH       = []; globalSolver = []; localSolver = []; 
   nTrial       = []; CLHMethod    = [];
else
   if isfield(Prob.CGO,'SCALE')
      SCALE = Prob.CGO.SCALE;
   else
      SCALE = [];
   end
   if isfield(Prob.CGO,'REPLACE')
      REPLACE = Prob.CGO.REPLACE;
   else
      REPLACE = [];
   end
   if isfield(Prob.CGO,'RandState')
      RandState = Prob.CGO.RandState;
   else
      RandState = [];
   end
   if isfield(Prob.CGO,'Percent')
      Percent = Prob.CGO.Percent;
   else
      Percent = [];
   end
   if isfield(Prob.CGO,'nSample')
      nSample = Prob.CGO.nSample;
   else
      nSample = [];
   end
   if isfield(Prob.CGO,'AddMP')
      AddMP = Prob.CGO.AddMP;
   else
      AddMP = [];
   end
   if isfield(Prob.CGO,'SMOOTH')
      SMOOTH = Prob.CGO.SMOOTH;
   else
      SMOOTH = [];
   end
   if isfield(Prob.CGO,'globalSolver')
      globalSolver = deblank(Prob.CGO.globalSolver);
   else
      globalSolver = [];
   end
   if isfield(Prob.CGO,'localSolver')
      localSolver = deblank(Prob.CGO.localSolver);
   else
      localSolver = [];
   end
   if isfield(Prob.CGO,'nTrial')
      nTrial = Prob.CGO.nTrial;
   else
      nTrial = [];
   end
   if isfield(Prob.CGO,'CLHMethod')
      CLHMethod = Prob.CGO.CLHMethod;
   else
      CLHMethod = [];
   end
end

if isempty(Prob.GO)
   GOMaxFunc = []; GOMaxIter     = []; maxFunc1 = [];  maxFunc2 = [];
   maxFunc3  = []; GOlocalSolver = []; DIRECT   = [];
else
   if isfield(Prob.GO,'MaxFunc')
      GOMaxFunc = Prob.GO.MaxFunc;
   else
      GOMaxFunc = [];
   end
   if isfield(Prob.GO,'MaxIter')
      GOMaxIter = Prob.GO.MaxIter;
   else
      GOMaxIter = [];
   end
   if isfield(Prob.GO,'maxFunc1')
      maxFunc1 = Prob.GO.maxFunc1;
   else
      maxFunc1 = [];
   end
   if isfield(Prob.GO,'maxFunc2')
      maxFunc2 = Prob.GO.maxFunc2;
   else
      maxFunc2 = [];
   end
   if isfield(Prob.GO,'maxFunc3')
      maxFunc3 = Prob.GO.maxFunc3;
   else
      maxFunc3 = [];
   end
   if isfield(Prob.GO,'localSolver')
      GOlocalSolver = Prob.GO.localSolver;
   else
      GOlocalSolver = [];
   end
   if isfield(Prob.GO,'DIRECT')
      DIRECT = Prob.GO.DIRECT;
   else
      DIRECT = [];
   end
end
% If not defined, assume problem is SMOOTH enough for local search
if isempty(SMOOTH), SMOOTH = 1; end

% Integer variables
IntVars  = DefPar(Prob.MIP,'IntVars',[]);

% Logical vector for integers
IV = false(d,1);

if isempty(IntVars)
   % No binary variables B or integer variables of type I
elseif any(IntVars==0) || all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > d)
      error('getCGOProbVars: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars = find(IV);
Reals   = find(~IV);

if SMOOTH
   if length(IntVars) == d
      %if isempty(globalSolver), globalSolver = 'glcDirect'; end
      if isempty(globalSolver), globalSolver = 'multiMINLP'; end
      % Set LOCAL below
   else
      if isempty(globalSolver), globalSolver = 'glcCluster'; end
      if strcmpi(globalSolver,'oqnlp') | strcmpi(globalSolver,'glcCluster') |...
         strcmpi(globalSolver,'multiMin') | strcmpi(globalSolver,'minlpBB')
         LOCAL = 0;
      else % Assume other solvers need a local search as well
         LOCAL = 1;
      end
   end
   if isempty(localSolver),  localSolver = GetSolver('con',1,0); end
else
   % What is the best non-smooth global solver???
   if isempty(globalSolver), globalSolver = 'glcDirect'; end
   LOCAL = 0; 
   % HKH Note - Here a robust non-smooth local solver should be used!!!
   if isempty(localSolver),  localSolver = GetSolver('con',0,0); end
end

if ~isempty(IntVars)
   % Safe guard bounds onto integer values
   x_L(IntVars)   = ceil(x_L(IntVars)); 
   x_U(IntVars)   = floor(x_U(IntVars)); 
   x_D            = x_U - x_L;
   if length(IntVars) == d
      % Pure IP problem
      nMax = prod(1+(x_U-x_L));
   else
      ix = true(d,1);
      ix(IntVars) = 0;
      if all(x_L(ix)==x_U(ix))
         % All continuous variables are fixed
         nMax = prod(1+(x_U(IntVars)-x_L(IntVars)));
      end
   end
   if ~isinf(nMax)
      % Pure IP problem
      LOCAL        = 0;
      switch lower(globalSolver)
        case {'multiminlp'}
           if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
           if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
        case {'glccluster'}
           globalSolver = glcDirect;
           fprintf('Pure IP problem, Set globalSolver = %s\n',globalSolver);
           if isempty(GOMaxFunc),    GOMaxFunc = min(nMax+1,100000); end
           if isempty(GOMaxIter),    GOMaxIter = 80000; end
        case {'glcdirect','glcfast','glcsolve'}
           if isempty(GOMaxFunc),    GOMaxFunc = min(nMax+1,100000); end
           if isempty(GOMaxIter),    GOMaxIter = 80000; end
        case {'glcsolve'}
           if isempty(GOMaxFunc),    GOMaxFunc = min(nMax+1,30000); end
           if isempty(GOMaxIter),    GOMaxIter = 20000; end
        case {'multimin'}
           if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
           if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
        case {'minlpbb'}
           if isempty(GOMaxFunc),    GOMaxFunc = max(10000,1000*d); end
           if isempty(GOMaxIter),    GOMaxIter = max(5000,500*d); end
        otherwise
           globalSolver = 'glcDirect';
           fprintf('Pure IP problem, Set globalSolver = %s\n',globalSolver);
           if isempty(GOMaxFunc),    GOMaxFunc = min(nMax+1,100000); end
           if isempty(GOMaxIter),    GOMaxIter = 80000; end
      end
   elseif strcmpi(globalSolver,'oqnlp')
      LOCAL = 0;
   elseif strcmpi(globalSolver,'glcDirect')
   elseif strcmpi(globalSolver,'glcSolve')
   elseif strcmpi(globalSolver,'minlpSolve')
   elseif strcmpi(globalSolver,'minlpBB')
      LOCAL = 0;
   elseif strcmpi(globalSolver,'glcCluster')
      LOCAL = 0;
   elseif strcmpi(globalSolver,'multiMin')
      LOCAL = 0;
   else % Default MINLP subsolver glcCluster
      globalSolver = 'glcCluster';
      LOCAL = 0;
   end
   SCALE = 0; % No scaling for MINLP or INLP-problems
end

% Default CGO input
if isempty(SCALE)
  if all(x_L==0) && all(x_U==1)
    SCALE = 0;
  else
    SCALE = 1;
  end
end
if isempty(RandState),    RandState = 0; end

if isempty(GOlocalSolver)
   GOlocalSolver  = localSolver;
end
if isempty(DIRECT)
   DIRECT  = 'glcDirect';
end

% MODIFICATION LOG:
%
% 080410 hkh Written 
% 080410 hkh Revise comments, set SMOOTH and DIRECT
% 080416 hkh Set LOCAL dependent on solver
% 080531 hkh Lower case solver names in Case statements, add MINLPBB for pureIP
