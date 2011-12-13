%function [Name, O, F, X, F_m, F00, Cc, nInit, Fpen, nFunc, n, MaxFunc] =...
%   CGOWarmStart(Prob, MaxFunc, PriLev)
%
% CGOWarmStart sets initial values on CGO variables if warm start,
% reading variables from Pro.CGO.WarmStartInfo, or from file cgoSave.mat
%

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written April 10, 2008.    Last modified August 24, 2009.

function [Name, O, F, X, F_m, F00, Cc, nInit, Fpen, nFunc, n, MaxFunc] =...
   CGOWarmStart(Prob, MaxFunc, PriLev)

% Restart with values from previous run.

WS = DefPar(Prob.CGO,'WarmStartInfo');
if ~isempty(WS)
   Name     = WS.Name;
   O        = WS.O;
   F        = WS.F;
   X        = WS.X;
   F_m      = WS.F_m;
   F00      = WS.F00;
   Cc       = WS.Cc;
   nInit    = WS.nInit;
   Fpen     = WS.Fpen;
   %fMinIdx  = WS.fMinIdx;
   rngState = WS.rngState;
else
   load('cgoSave.mat','Name','O','F','X','F_m','F00','Cc','nInit','Fpen',...
        'rngState');
        %'fMinIdx','rngState');
end
rand('state',rngState);

Name1 = deblank(Prob.Name);   % Name for the problem

if strcmp(Name1,Name)
   [d n] = size(X);
   nFunc = n; % Total count of function evaluations

   if PriLev > 0
      fprintf('\n Restarting with %d sampled points from previous run\n',n);
   end
else
   nFunc = 0; % Total count of function evaluations
   n     = 0; 
   if PriLev >= -1000
      fprintf('Previous run was with Problem %s\n',Name);
      fprintf('This run is with Problem %s\n',Name1);
      fprintf('Impossible to do restart.\n');
      fprintf('Maybe there exists several files cgoSave.mat?\n');
   end
end

% Add extra function evaluations 
if MaxFunc <= nFunc
   MaxFunc = MaxFunc + nFunc;
end
   
% MODIFICATION LOG:
%
% 080410 hkh Written 
% 080414 hkh Revised comments 
% 080701 hkh Added F00 and Cc, to handle costly constraints
% 090824 hkh No need to define fMinIdx - not returned
