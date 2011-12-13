% function [f,Result] = sepMIPD_f(x, Prob)
%
% Function computing NLP subproblem, when using standard TOMLAB routines,
% e.g. glcSolve and glcFast for IP part
%
% Global variables used: counter fBest
%
% sepMIPD_f is setting up the shrinked NLP problem.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Apr 10, 2004. Last modified July 25, 2011.

function [f,Result] = sepMIPD_f(x, Prob)

global counter fBest
counter          = counter + 1;
P                = Prob.Prob;
IntVars          = P.MIP.IntVars;
BLENDS           = P.user.BLENDS;
RM               = P.user.RM;
ix               = find(x);
% Create a smaller problem
n                = length(ix);
N                = n*BLENDS;
ixF = [];
v   = RM;
for i = 1:BLENDS
    ixF = [ixF;v+ix];
    v = v + RM;
end
P.user.c         = Prob.user.c(ixF);
if ~isinf(fBest)
   c1 = sum(P.user.c*0.03);
   c2 = min(P.user.c)*(1-0.03*n);
   fNew = c1 + c2;
   if fNew > fBest
      fprintf('SKIP Point %d, best f(x) possible: %f. ',fNew);
      fprintf('Current best f(x):  %f\n',fNew);
      f      = fNew;
      Result = [];
      return
   end
end
P.N              = N;
P.x_L            = 0.03*ones(N,1);
P.x_U            = ones(N,1);
z                = 1/n;
P.x_0            = z*ones(N,1);
P.ConsPattern    = P.ConsPattern(:,ixF);
if ~isempty(P.LINCON)
   % Additional linear constraints added last for nonlinear part
   P.A     = [P.A  ;P.LINCON.A(:,ixF)];
   P.b_L   = [P.b_L;P.LINCON.b_L];
   P.b_U   = [P.b_U;P.LINCON.b_U];
   P.mLin  = size(P.A,1);
end
P.user.blx       = Prob.user.blx(ix);
P.user.tgNx      = Prob.user.tgNx(ix);
P.user.tgGx      = Prob.user.tgGx(ix);
P.user.safax     = Prob.user.safax(ix);
P.user.pufax     = Prob.user.pufax(ix);
P.user.h2mx      = Prob.user.h2mx(ix);
P.user.Nmodels   = Prob.user.Nmodels;
P.user.RM        = n;
P.user.RMIP      = 0;
P.ConIx          = findpatt(P.ConsPattern);

P.P   = counter;
P.P2  = counter;
P.Pnr = Prob.P;

[f,Result] = nlpblend(ix, P);

h_k = Result.h_k;
if h_k < 1E-3 & f < fBest
   fBest = f;
end
f               = f + 1000000*h_k;
Result.IntVars = ix;

% MODIFICATION LOG:
%
% 070519  hkh  Written
% 110725  hkh  ConsIx changed to ConIx
