% mmx_dc.m
%
% function dc=mmx_dc(x, Prob, varargin)
%
% mmx_dc computes the gradient of the constraints c at the point x
% and the Jacobian matrix for the minimax residuals r(x).
% r is the residuals in the original formulation: min max r(x)
% The extra variable x(n+1) has derivatives -1 exactly for each residual.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written May 22, 1999. Last modified July 22, 2011.

function dc=mmx_dc(x, Prob, varargin)

global n_dc NLP_xc NLP_c NLP_xdc NLP_dc mad_c mad_dc
global mad_r
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

J             = nlp_J(x(1:n+1),Prob, varargin{:});

if isempty(Prob.FUNCS.c)
   dc = [];
else
   % Must save and reset global variables
   n_dc1         = n_dc;
   NARG1         = NARG;
   NLP_xc1       = NLP_xc;
   NLP_c1        = NLP_c;
   n_dc          = 0;
   NARG          = args;
   NLP_xdc1      = NLP_xdc;
   NLP_dc1       = NLP_dc;
   mad_c1        = mad_c;
   mad_dc1       = mad_dc;

   Prob.ConsDiff = Prob.minimax.ConsDiff;
   Prob.ADCons   = Prob.minimax.ADCons;
   % Pick up these 3 global variables from mmx_c call
   global NLP_xcS NLP_cS mad_cS
   NLP_xc        = NLP_xcS;
   NLP_c         = NLP_cS;
   mad_c         = mad_cS;

   NLP_xdc       = [];
   NLP_dc        = [];
   dc            = nlp_dc(x(1:n+1),Prob, varargin{:});
   n_dc          = n_dc1;
   NARG          = NARG1;
   NLP_xc        = NLP_xc1;
   NLP_c         = NLP_c1;
   NLP_xdc       = NLP_xdc1;
   NLP_dc        = NLP_dc1;
   mad_c         = mad_c1;
   mad_dc        = mad_dc1;
end

if isempty(dc)
   dc=sparse([J,-ones(size(J,1),1)]);
else
   dc=sparse([dc,zeros(size(dc,1),1);[J,-ones(size(J,1),1)]]);
end

% MODIFICATION LOG
%
% 990522  hkh  minimax constraint derivative routine written
% 020328  hkh  Modified for Tomlab v3.1
% 020404  hkh  Calling all nlp_XXX routines with full x vector
% 020411  hkh  Reset global variables before calling nlp_dc
% 020416  hkh  Set original Prob.NumDiff and Prob.ConsDiff before call
% 040126  hkh  Wrong global variable, n_dc not n_c should be used
% 041130  hkh  Wrong check for nonlinear constraints, check Prob.FUNCS.c not dc
% 060814  med  FUNCS used for callbacks instead
% 110722  hkh  Saving 3 globals in mmx_c to handle AD and ConsDiff, used here 
