%function [X,O,F,Fpen,F00,Cc,C,nCon,nSample,ExDText,initRS] = expDesignD(...
%   Percent, nSample, AddMP, nTrial, CLHMethod, SCALE, RandState, ...
%   PriLev, Prob, varargin)
%
% expDesignD calls expDesign, after safe guarding some input that may
% not be set OK if called directly from command line
%
% expDesign computes an initial experimental design for CGO solvers.
%
% expDesignD has the same input and output as expDesign
% See the help for expDesign
%

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., $Release: 7.4.0$
% Written April 10, 2008.    Last modified August 24, 2008.


function [X,O,F,Fpen,F00,Cc,C,nCon,nSample,ExDText,initRS] = expDesignD(...
   Percent, nSample, AddMP, nTrial, CLHMethod, SCALE, RandState, ...
   PriLev, Prob, varargin)

if nargin < 9
   error('expDesignD needs 9 parameters');
end

% Safe guard a bit, if expDesign called directly
if isempty(Prob.optParam)
   bTol = 1E-7;
   cTol = 1E-5;
else
   if isfield(Prob.optParam,'bTol')
      bTol = Prob.optParam.bTol;      % Linear constraint feasibility tolerance
      % if isempty(bTol), bTol = 1E-7; end
   else
      bTol = 1E-7;
   end
   if isfield(Prob.optParam,'cTol')
      cTol = Prob.optParam.cTol;      % Constraint feasibility tolerance
      % if isempty(cTol), cTol = 1E-5; end
   else
      cTol = 1E-5;
   end
end
Prob.optParam.bTol = bTol; 
Prob.optParam.cTol = cTol;

%d= length(x_D);

dLin = size(Prob.A,1);
if dLin > 0
   Prob.b_L = [Prob.b_L;-inf*ones(dLin-length(Prob.b_L),1)];
   Prob.b_U = [Prob.b_U; inf*ones(dLin-length(Prob.b_U),1)];
end
dCon = max(length(Prob.c_L),length(Prob.c_U));
if dCon > 0
   Prob.c_L = [Prob.c_L;-inf*ones(dCon-length(Prob.c_L),1)];
   Prob.c_U = [Prob.c_U; inf*ones(dCon-length(Prob.c_U),1)];
end

if 0
%NHQ   Fix when called standalone
%
% Check boundries. If one empty, set to +/-inf
flag(1) = isempty(Prob.b_L);
flag(2) = isempty(Prob.b_U);
flag(3) = isempty(Prob.c_L);
flag(4) = isempty(Prob.c_U);

if sum(flag(1:2)) == 1
   if flag(1) == 1
      Prob.b_L = -inf;
   else
      Prob.b_U =  inf;
   end
end
if sum(flag(3:4)) == 1
   if flag(3) == 1
      Prob.c_L = -inf;
   else
      Prob.c_U =  inf;
   end
end
end

if ~isfield(Prob,'NARG')
   NARG=zeros(11,1);
   z= Prob.FUNCS.f;
   if ~isempty(z), NARG(1)  = xnargin(z); end
   z= Prob.FUNCS.c;
   if ~isempty(z), NARG(4)  = xnargin(z); end
   z= Prob.FUNCS.dc;
   if ~isempty(z), NARG(5)  = xnargin(z); end
   z= Prob.FUNCS.d2c;
   if ~isempty(z), NARG(6)  = xnargin(z); end
   z= Prob.FUNCS.fc;
   if ~isempty(z), NARG(10) = xnargin(z); end
   Prob.NARG = NARG;
end

[X,O,F,Fpen,F00,Cc,C,nCon,nSample,ExDText,initRS] = expDesign(Percent, ...
   nSample,AddMP, nTrial, CLHMethod, SCALE, RandState, PriLev, Prob, varargin);

% MODIFICATION LOG:
%
% 080410 hkh Written 
% 080417 hkh Add output variable ExDText
% 080701 hkh Modified because of changes in expDesign
% 080717 hkh Added definition of some elements in Prob.NARG, if not defined
% 090724 hkh Must set Prob.NARG = NARG
