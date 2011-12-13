% L2_dc.m
%
% function dc=L2_dc(x, Prob, varargin)
%
% L2_dc computes the gradient of the constraints c at the point x
% and the Jacobian matrix for the L2 residuals r(x).
% r is the residuals in the original formulation: min 0.5*r'*r
% The extra variables have derivatives -1 exactly


% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Apr 13, 2002. Last modified July 22, 2011.

function dc=L2_dc(x, Prob, varargin)

global n_dc NLP_xc NLP_c NLP_xdc NLP_dc mad_c mad_dc
global mad_r
args          = Prob.L2.args;
global NARG
NARG(7:8)     = args(7:8);

m             = Prob.L2.m;
n             = length(x) - m;

Prob.x_0      = Prob.x_0(1:n); % Adjust x_0 for user routine
Prob.x_L      = Prob.x_L(1:n);
Prob.x_U      = Prob.x_U(1:n);
Prob.N        = n;

Prob.FUNCS.r  = Prob.L2.r;
Prob.FUNCS.J  = Prob.L2.J;
Prob.FUNCS.c  = Prob.L2.c;
Prob.FUNCS.dc = Prob.L2.dc;
Prob.NumDiff  = Prob.L2.NumDiff;
Prob.ADObj    = Prob.L2.ADObj;

J             = nlp_J(x,Prob, varargin{:});

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

   Prob.ConsDiff = Prob.L2.ConsDiff;
   Prob.ADCons   = Prob.L2.ADCons;
   % Pick up these 3 global variables from L1_c call
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
   dc=sparse([J,-speye(m,m)]);
else
   dc=sparse([dc,zeros(size(dc,1),m);[J,-speye(m,m)]]);
end

% MODIFICATION LOG
%
% 020413  hkh  Written
% 020416  hkh  Set original Prob.NumDiff and Prob.ConsDiff before call
% 040126  hkh  Wrong global variable, n_dc not n_c should be used
% 040126  hkh  Field Prob.L1.ConsDiff should be Prob.L2.ConsDiff
% 041130  hkh  Wrong check for nonlinear constraints, check Prob.FUNCS.c not dc
% 060814  med  FUNCS used for callbacks instead
% 110722  hkh  Saving 3 globals in L1_c to handle AD and ConsDiff, used here 
