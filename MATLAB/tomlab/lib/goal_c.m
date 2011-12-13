% goal_c.m
%
% function cx=goal_c(x, Prob, varargin)
%
% goal_c computes the goal attainment constraints c and r in the point x
% r is the residual functions in the original formulation

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Sept 28, 2002. Last modified Aug 14, 2006.

function cx=goal_c(x, Prob, varargin)

global n_c NLP_xc NLP_c

n=length(x)-1;

Prob.x_0 = Prob.x_0(1:n); % Adjust x_0 for user routine
Prob.x_L = Prob.x_L(1:n);
Prob.x_U = Prob.x_U(1:n);
Prob.N   = n;

args          = Prob.minimax.args;
Prob.FUNCS.r  = Prob.minimax.r;
Prob.FUNCS.J  = Prob.minimax.J;
Prob.ConsDiff = Prob.minimax.ConsDiff;
Prob.NumDiff  = Prob.minimax.NumDiff;

global NARG
NARG(7:8) = args(7:8);

% Send also x(n+1) to nlp_r. nlp_r calls user routine only with x(1:n),
% where n is Prob.N

r  = nlp_r(x(1:n+1),Prob, varargin{:});

Prob.FUNCS.c  = Prob.minimax.c;
Prob.FUNCS.dc = Prob.minimax.dc;

if isempty(Prob.FUNCS.c)
    cx = [];
else
    % Must save and reset global variables
    n_c1    = n_c;
    NARG1   = NARG;
    NLP_xc1 = NLP_xc;
    NLP_c1  = NLP_c;
    n_c     = 0;
    NARG    = args;
    NLP_xc  = [];
    NLP_c   = [];
    cx = nlp_c(x(1:n+1),Prob, varargin{:});

    n_c     = n_c1;
    NARG    = NARG1;
    NLP_xc  = NLP_xc1;
    NLP_c   = NLP_c1;
end

cx = [cx;r-x(n+1)*Prob.LS.w];

% MODIFICATION LOG
%
% 020928  hkh  Created based on mmx_c.m
% 030127  hkh  Must set Prob.ConsDiff = Prob.minimax.ConsDiff and
% 030127  hkh  Prob.NumDiff  = Prob.minimax.NumDiff (if opt tbx interface)
% 040126  hkh  Missing semi colon
% 060814  med  FUNCS used for callbacks instead