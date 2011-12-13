% mmx_c.m
%
% function cx=mmx_c(x, Prob, varargin)
%
% mmx_c computes the minimax constraints c and r in the point x
% r is the residuals in the original formulation: min max r(x)

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written May 22, 1999. Last modified July 22, 2011.

function cx=mmx_c(x, Prob, varargin)

global n_c NLP_xc NLP_c mad_c

% Globals for nlp_r and nlp_J need not be reset or saved
% They are unique and not set by mmx_f, mmx_g and mmx_H
% global mad_r is used in mmx_dc to compute AD Jacobian
global wLS n_J LS_x LS_r LS_xJ LS_J mad_r mad_J

args          = Prob.minimax.args;
global NARG
NARG(7:8)     = args(7:8);

n             = length(x)-1;
Prob.x_0      = Prob.x_0(1:n); % Adjust x_0 for user routine
Prob.x_L      = Prob.x_L(1:n);
Prob.x_U      = Prob.x_U(1:n);
Prob.N        = n;

Prob.FUNCS.r  = Prob.minimax.r;
Prob.FUNCS.J  = Prob.minimax.J;
Prob.FUNCS.c  = Prob.minimax.c;
Prob.FUNCS.dc = Prob.minimax.dc;
Prob.NumDiff  = Prob.minimax.NumDiff;
Prob.ADObj    = Prob.minimax.ADObj;

% Send also x(n+1) to nlp_r. nlp_r calls user routine only with x(1:n),
% where n is Prob.N

r             = nlp_r(x(1:n+1),Prob, varargin{:});

if isempty(Prob.FUNCS.c)
   cx = [];
else
   % Must save and reset global variables
   % These 3 globals are then used in mmx_dc
   global NLP_xcS NLP_cS mad_cS

   n_c1          = n_c;
   NARG1         = NARG;
   NLP_xc1       = NLP_xc;
   NLP_c1        = NLP_c;
   mad_c1        = mad_c;

   n_c           = 0;
   NARG          = args;
   NLP_xc        = [];
   NLP_c         = [];
   Prob.ConsDiff = Prob.minimax.ConsDiff;
   Prob.ADCons   = Prob.minimax.ADCons;
   cx            = nlp_c(x(1:n+1),Prob, varargin{:});
   % Save these three globals
   NLP_xcS       = NLP_xc;
   NLP_cS        = NLP_c;
   mad_cS        = mad_c;
   % Reset these 5 globals
   n_c           = n_c1;
   NARG          = NARG1;
   NLP_xc        = NLP_xc1;
   NLP_c         = NLP_c1;
   mad_c         = mad_c1;
end

cx = [cx;r-x(n+1)];

% MODIFICATION LOG
%
% 990522  hkh  minimax constraint routine written
% 020404  hkh  Calling all nlp_XXX routines with full x vector
% 020411  hkh  Save global before calling nlp_c
% 040126  hkh  Missing semi colon
% 060814  med  FUNCS used for callbacks instead
% 110722  hkh  Saving 3 globals to handle AD and ConsDiff, use in mmx_dc
