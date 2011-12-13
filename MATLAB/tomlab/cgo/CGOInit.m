%function [Name, O, F, X, F_m, F00, Cc, nInit, Fpen, nFunc, n, nCon,nSample,...
%   ExDText] = CGOInit(Prob, RandState, nSample, Percent, AddMP, nTrial,...
%   CLHMethod, SCALE, REPLACE, PriLev, varargin)
%
% CGOInit sets initial values on CGO variables if no warm start
%

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written April 10, 2008.    Last modified July 1, 2008.

function [Name, O, F, X, F_m, F00, Cc, nInit, Fpen, nFunc, n, nCon,nSample,...
   ExDText] = CGOInit(Prob, RandState, nSample, Percent, AddMP, nTrial,...
   CLHMethod, SCALE, REPLACE, PriLev, varargin)

Name = deblank(Prob.Name);  % Problem name

% Set pseudo random generator
if isnan (RandState)
   % Do no initialization
elseif length(RandState) >1 || RandState >= 0
   rand('state',RandState);
else
   rand('state',sum(100*clock));
end

% Send NaN as RandState to expDesign, in order to not change its value.
% Many DIFFERENT INITIALIZATION STRATEGIES, plus user defined
if Prob.simType == 1
   % Constraints are costly, avoid use in expDesign
   fprintf('Constraints are costly, avoid use in expDesign\n');
elseif Prob.simType == 2
   fprintf('Constraints both noncostly (used in expDesign) and costly\n');
end

[X,O,F,Fpen,F00,Cc,C,nCon,nSample,ExDText,initRS] = expDesign(Percent, ...
   nSample, AddMP, nTrial, CLHMethod, SCALE, NaN, PriLev, Prob, varargin{:});

% HKH Note - initRS is never used

% Set dimension d and number of points n
[d n] = size(X);
  
% Set initial number of points nInit to n
nInit = n;
nFunc = n; % Number of function evaluations.

if REPLACE > 1
   % Replace very large function values by log
   FMIN = min(F);
   if FMIN <= 0
      FMAX = 10^REPLACE;
   else
      FMAX = 10^(ceil(log10(FMIN))+REPLACE);
   end
   ix = find(F>FMAX);
   F_m = F;
   if ~isempty(ix)
      fprintf('REPLACE %d large values with ',length(ix));
      fprintf('FMAX+log10(F-FMAX+1). FMAX=%e\n',FMAX);
      F_m(ix) = FMAX + log10(F(ix)-FMAX+1);
   end
elseif REPLACE == 1
   % Replace large function values by median(F)
   F_m = min(median(F),F);
else
   F_m = F;
end


% MODIFICATION LOG:
%
% 080410 hkh   Written 
% 080417 hkh   Add additional output ExDText
% 080617  frhe Made REPLACE>1 transformation monotonic
% 080617  frhe log10 used in REPLACE>1 in accordance with help
% 080630  hkh  Handle both costly and noncostly constraints. Add Cc as output
% 080701  hkh  Add F00, F without penalty from costly Cc
