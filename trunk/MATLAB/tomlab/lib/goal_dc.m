% goal_dc.m
%
% function dc=goal_dc(x, Prob, varargin)
%
% goal_dc computes the gradient of the constraints c at the point x
% and the Jacobian matrix for the goal attainment residual functions r(x).
% The extra variable x(n+1) has derivatives -w 

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Sept 28, 2002. Last modified Aug 14, 2006.

function dc=goal_dc(x, Prob, varargin)

global n_dc NLP_xc NLP_c

n=length(x)-1;
Prob.x_0 = Prob.x_0(1:n); % Adjust x_0 for user routine
Prob.x_L = Prob.x_L(1:n);
Prob.x_U = Prob.x_U(1:n);
Prob.N   = n;
args     = Prob.minimax.args;

global NARG
NARG(7:8) = args(7:8);

Prob.FUNCS.r = Prob.minimax.r;
Prob.FUNCS.J = Prob.minimax.J;

Prob.NumDiff = Prob.minimax.NumDiff;
J = nlp_J(x(1:n+1),Prob, varargin{:});

Prob.FUNCS.c  = Prob.minimax.c;
Prob.FUNCS.dc = Prob.minimax.dc;

if isempty(Prob.FUNCS.c)
   dc = []; % No nonlinear constraints
else
   % Must save and reset global variables
   n_dc1   = n_dc;
   NARG1   = NARG;
   NLP_xc1 = NLP_xc;
   NLP_c1  = NLP_c;
   n_dc    = 0;
   NARG    = args;
   NLP_xc  = [];
   NLP_c   = [];
   Prob.ConsDiff = Prob.minimax.ConsDiff;

   dc = nlp_dc(x(1:n+1),Prob, varargin{:});
   n_dc   = n_dc1;
   NARG   = NARG1;
   NLP_xc = NLP_xc1;
   NLP_c  = NLP_c1;
end

if isempty(dc)
   dc=sparse([J,-Prob.LS.w]);
else
   dc=sparse([dc,zeros(size(dc,1),1);[J,-Prob.LS.w]]);
end

% MODIFICATION LOG
%
% 020928  hkh  Created based on mmx_c.m
% 030611  ango Changed Prob.user.w (and .g) to Prob.LS.w (and g),
%              to match information in goalSolve.m 
% 040126  hkh  Wrong global variable, n_dc not n_c should be used
% 041124  hkh  Wrong check for nonlinear constraints, check Prob.FUNCS.c not dc
% 060814  med  FUNCS used for callbacks instead