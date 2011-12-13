% L1_c.m
%
% function cx=L1_c(x, Prob, varargin)
%
% L1_c computes the L1 constraints c and r in the point x
% r is the residuals in the original formulation: min sum |r_i(x)|

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Apr 13, 2002. Last modified July 22, 2011.

function cx=L1_c(x, Prob, varargin)

global n_c NLP_xc NLP_c mad_c

% Globals for nlp_r and nlp_J need not be reset or saved
% They are unique and not set by L1_f, L1_g and L1_H
% global mad_r is used in L1_dc to compute AD Jacobian
global wLS n_J LS_x LS_r LS_xJ LS_J mad_r mad_J

args          = Prob.L1.args;
global NARG
NARG(7:8)     = args(7:8);

m             = Prob.L1.m;
n             = length(x) - 2*m;

Prob.FUNCS.r  = Prob.L1.r;
Prob.FUNCS.J  = Prob.L1.J;
Prob.FUNCS.c  = Prob.L1.c;
Prob.FUNCS.dc = Prob.L1.dc;
Prob.NumDiff  = Prob.L1.NumDiff;
Prob.ADObj    = Prob.L1.ADObj;

Prob.x_0      = Prob.x_0(1:n); % Adjust x_0 for user routine
Prob.x_L      = Prob.x_L(1:n);
Prob.x_U      = Prob.x_U(1:n);
Prob.N        = n;

% Send also x(n+1:n+2*m) to nlp_r. nlp_r calls user routine only with x(1:n),
% where n is Prob.N

r             = nlp_r(x,Prob, varargin{:});
m             = length(r);

if isempty(Prob.FUNCS.c)
   cx = [];
else
   % Must save and reset global variables
   % These 3 globals are then used in L1_dc
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
   Prob.ConsDiff = Prob.L1.ConsDiff;
   Prob.ADCons   = Prob.L1.ADCons;
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

cx = [cx;r+x(n+1:n+m)-x(n+m+1:n+2*m)];

% MODIFICATION LOG
%
% 020413  hkh  Written
% 040126  hkh  Missing semi colon
% 060814  med  FUNCS used for callbacks instead
% 110722  hkh  Saving 3 globals to handle AD and ConsDiff, use in L1_dc
