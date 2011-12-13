% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written Nov 24, 2006. Last modified Nov 24, 2006.

function gnG = gnGrid(Result)

if nargin < 1
   error('gnGrid needs one parameter, the structure Result');
end

epsX      = 1.7E-6;
PriLev    = 0;

snProb = Result.CGO.snProb;
gnProb = Result.CGO.gnProb;

% Pick up all parameters
Name    = Result.CGO.WarmStartInfo.Name;
nInit   = Result.CGO.WarmStartInfo.nInit;
fMinIdx = Result.CGO.WarmStartInfo.fMinIdx;


RandState    = Result.CGO.RandState;
Percent      = Result.CGO.Percent;
nSample      = Result.CGO.nSample;
fTol         = Result.CGO.fTol;
infStep      = Result.CGO.infStep;
N            = Result.CGO.N;
AddMP        = Result.CGO.AddMP;
idea         = Result.CGO.idea;
rbfType      = Result.CGO.rbfType;
fStarRule    = Result.CGO.fStarRule;
eps_sn       = Result.CGO.eps_sn;
DeltaRule    = Result.CGO.DeltaRule;
REPLACE      = Result.CGO.REPLACE;
SCALE        = Result.CGO.SCALE;
MaxCycle     = Result.CGO.MaxCycle;
LOCAL        = Result.CGO.LOCAL;
localSolver  = Result.CGO.localSolver;
globalSolver = Result.CGO.globalSolver;
globalSolver = 'oqnlp';
localSolver  = 'oqnlp';
GOMaxFunc    = Result.CGO.GOMaxFunc;
GOMaxIter    = Result.CGO.GOMaxIter;
maxFunc1     = Result.CGO.maxFunc1;
maxFunc2     = Result.CGO.maxFunc2;
maxFunc3     = Result.CGO.maxFunc3;
fGoal        = Result.CGO.fGoal;

x_L          = Result.Prob.x_L;
x_U          = Result.Prob.x_U;
x_D          = x_U - x_L;
d            = length(x_D);
IntVars      = Result.Prob.MIP.IntVars;
IV           = false(d,1);
IV(IntVars)  = 1;
Reals        = find(~IV);

% Use optimal value for printing, if given
if isempty(snProb.x_opt)
   x_opt = [];
else
   x_opt = snProb.x_opt(1,:);
end
if SCALE
   if ~isempty(x_opt)
      xOptS   = (x_opt(:)-x_L)./x_D; 
   else
      xOptS   = [];
   end
else
   xOptS   = x_opt(:);
end

DEBUG        = 0;
Its          = Result.CGO.Its;

nPnt = length(Result.CGO.WarmStartInfo.F);

nPnt = min(15,nPnt);
nPnt = (16:25);

for jj=1:length(nPnt)

n = nPnt(jj);
fprintf('\nTarget value grid using %d interpolation points. ',n);
fprintf('Initial set %d points',nInit);
fprintf('\n');

O       = Result.CGO.WarmStartInfo.O(:,1:n);
F       = Result.CGO.WarmStartInfo.F(1:n);
X       = Result.CGO.WarmStartInfo.X(:,1:n);
F_m     = Result.CGO.WarmStartInfo.F_m(1:n);
Fpen    = Result.CGO.WarmStartInfo.Fpen(1:n);

MaxFunc = size(X,2) + 1;

% TOMSOL INIT, send F_m to Fortran, not F
control = tomsol(27, MaxFunc, X, F_m, rbfType, idea, DEBUG, REPLACE);
if control < 0
    fprintf('Initial interpolation failed');
    tomsol(25); % Deallocates memory
    return      %Something is really wrong
end

z = Fpen;
% Set infeasible points as infinity before check on minimum
z(Fpen-F >= 1E-14)=Inf;
   
[fMin,fIdx] = min(z);

if isinf(fMin)
   % No feasible point found
   Feasible = 0;
   % Take the best infeasible point
   [fMin,fIdx] = min(Fpen);
else
   Feasible = 1;
end

% Best point found in unit space
x_min = X(:,fIdx(1)); 


% Best point found in original space
if SCALE
   O_min = tomsol(9, x_L, x_min, x_D); 
else
   O_min = x_min;
end

% Set parameters used in snSolve
snProb.PriLev           = PriLev;
snProb.epsX             = epsX;
if ~strcmpi(globalSolver,localSolver)
   snProb.dLin             = dLin;
   snProb.MaxIter          = MaxIter;
end

snProb.CGO        = gnProb.CGO;
% Set parameters used in snSolve
snProb.CGO.globalSolver = globalSolver;
snProb.CGO.localSolver  = localSolver;
snProb.CGO.fnStar = NaN;
snResult          = snSolve(snProb);
min_sn            = snResult.f_k;
if isempty(snResult.x_k)
   min_sn_y = x_min;
else
   min_sn_y = snResult.x_k(:,1);
end

% Set parameters used in gnSolve
% Set parameters in global structure CGO
gnProb.CGO.X            = X;
gnProb.CGO.globalSolver = globalSolver;
gnProb.PriLev           = PriLev;
gnProb.PriLevOpt        = 0;
gnProb.LOCAL            = LOCAL;
gnProb.epsX             = epsX;
if LOCAL
   gnProb.CGO.localSolver  = localSolver;
   gnProb.dLin             = dLin;
   gnProb.MaxIter          = MaxIter;
end

% Set iteration dependent parameters used in gnSolve, output from snSolve
% gnProb.IX   = snResult.multiMin.IX;
% gnProb.xOpt = snResult.multiMin.xOpt;
gnProb.IX   = [];
gnProb.xOpt = [];

modN = -9;
gnProb.CGO.modN   = modN;
v = [1E-8 1E-6 1E-4 1E-3 1E-2 1E-1 0.25 0.5 1 1.5 2 4 8 12 16 20 100 1000];
fnS = min_sn-v*abs(min_sn);
xprint(fnS,'fnStar:');
xprint(fMin,'fMin:  ');

for i=1:length(fnS)
    fnStar = fnS(i);
    gnProb.CGO.fnStar = fnStar;

    gnResult = gnSolve(gnProb);

    if strcmpi(globalSolver,'multiMin')
       xNew   = gnResult.multiMin.xOpt(:,1);
       fNew   = gnResult.multiMin.fOpt(1);
       onBx   = gnResult.multiMin.onB;
       ix     = find(onBx==0);
       if ix(1) > 1
          doX    = min(tomsol(30,xNew,gnResult.multiMin.xOpt(:,ix(1))));
          fprintf('fOpt %11.6f, #%d, ',fNew,ix(1))
          fprintf('fInterior %f, ',gnResult.multiMin.fOpt(ix(1)));
          fprintf('Dist: xNew-xInt %f',doX);
          fprintf('\n');
       end
    else
       xNew   = gnResult.x_k(:,1);
       fNew   = gnResult.f_k;
    end
    Update = 0; % Not updated X yet
    % onB  = nOnBound(x_k(:,1),gnProb.x_L,gnProb.x_U,epsX);
    if i == 1
       [onB_sn, doX_sn, doM_sn] = statGN(snProb.x_L,snProb.x_U,min_sn_y,...
           x_min,X,epsX,Update,IntVars,Reals);
       if ~isempty(x_opt)
          doO_sn = min(tomsol(30,min_sn_y,xOptS));
       else
          doO_sn = [];
       end
       fprintf('min_sn %9.4f ',min_sn);
       fprintf('[%d] ',onB_sn);
       fprintf('doX %f ',doX_sn);
       fprintf('doM %f ',doM_sn);
       if ~isempty(x_opt)
          fprintf('doO %f ',doO_sn);
       end
       fprintf('\n');
    end
    % Distance between new point and best point found or closest point in X
    [onB, doX, doM] = statGN(gnProb.x_L,gnProb.x_U,xNew,x_min,X,epsX,...
       Update,IntVars,Reals);
    if ~isempty(x_opt)
       doO = min(tomsol(30,xNew,xOptS));
    else
       doO = [];
    end
    % Distance between new point and min on RBF surface
    SoO = min(tomsol(30,xNew,min_sn_y));

    fprintf('fnStar %9.4f ',fnStar);
    fprintf('fNew %11.6f ',fNew);
    fprintf('xNew: [%d] ',onB);
    fprintf('doX %f ',doX);
    fprintf('doM %f ',doM);
    if ~isempty(x_opt)
       fprintf('doO %f ',doO);
    end
    fprintf('SoO %f ',SoO);
    % New point in original space
    if SCALE
       O_new   = tomsol(9, x_L, xNew, x_D); 
    else
       O_new   = xNew;
    end
   % global NLP_x NLP_f NARG
   NLP_x=[]; NLP_f=[]; NARG = [];
   % fx  = nlp_f(O_new, Result.Prob, varargin{:});
   fx  = nlp_f(O_new, Result.Prob);
   fprintf('f(xNew) %f ',fx);
   fprintf('\n');
end

gnG = gnResult;

tomsol(25) % Deallocates memory

end

