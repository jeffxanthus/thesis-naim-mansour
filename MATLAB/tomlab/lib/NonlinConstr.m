% NonlinConstr computes constraint information using Prob structure info
%
% function [m, c_k, dc_k, cEqual, c_L, c_U] = NonlinConstr(Prob, varargin)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., Sweden. $Release: 7.3.0$
% Written June 28, 1999. Last modified Aug 13, 2009.

function [m, c_k, dc_k, cEqual, c_L, c_U] = NonlinConstr(Prob, varargin)

if nargin < 1
   error('NonlinConstr must have 1 parameter Prob as input');
end

% Nonlinear constraints

c_L = Prob.c_L(:);
c_U = Prob.c_U(:);
x_k = Prob.x_0;  

if isempty(c_L) && isempty(c_U)
   c_k=[];
else
   c_k=nlp_c(x_k, Prob, varargin{:});  % Constraints 
   c_k=c_k(:);
end

m = length(c_k);                  % # of constraints
n = length(x_k);                  % # of variables

if m == 0
   c_k   = zeros(0,1); % Because c_L(:),c_U(:) makes empty vector size (0,1)
   dc_k  = zeros(0,n); 
   cEqual= zeros(0,1);
else
   dc_k  = nlp_dc(x_k, Prob, varargin{:});

   if isempty(c_U),c_U= Inf*ones(m,1); end
   if isempty(c_L),c_L=-Inf*ones(m,1); end

   cEqual=eq(c_L,c_U);
end

% MODIFICATION LOG:
%
% 040506  hkh  Safeguard c_k, changing to a column vector
% 090813  med  mlint check