% gkls_prob: Defines global optimization problems with simple bounds.
%
% function [probList, Prob] = gkls_prob(P, ask, Prob);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure
%
% The TOMLAB GKLS implementation is based on the GKLS problem generator 
% by M. GAVIANO, D.E. KVASOV, D. LERA, and Ya.D. SERGEYEV.
% http://si.deis.unical.it/~yaro/GKLS.html

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.2.0$
% Written Sep 1, 1998.   Last modified Jun 5, 2008.

function [probList, Prob] = gkls_prob(P, varargin)

if nargin < 1
    P=[];
end

probList = [];

if isempty(P)
   return
end

% LATER: add the possibility for the ask parameter to control the following parameters.
% There's an infinite number of problems available, really

n        = 2;     % problem dimension, must be at least 2
nmin     = 4;     % number of minimizers wanted, at least 2
f_star   = -2.0;  % desired objective value at _global_ minimum
r_star   = 0.01;  % distance between global optimum and paraboloid vertex
rho_star = 0.005; % global minimum attraction region radius

nf = P;

x_L = -1.0*ones(n,1);
x_U =  1.0*ones(n,1);

% Tgkls checking is not complete:
% Some additional rules, which are NOT checked inside Tgkls. This way isn't
% very nice to do it, better to check this in some m-file and stop if not 
% satisfied:

% 1st requirement: 0 < r_star <  0.5*min|xu-xl|
d_min = min(x_U-x_L);

while r_star >= 0.5*d_min
   r_star = 0.75 * r_star;
end

% .. and furthermore, it is necessary that 0 < rho_star < 0.5*r_star
while rho_star >= 0.5*r_star
   rho_star = 0.75 * rho_star;
end

if P <= 100
   % ND problems, only f(x) is available
   f = 'gkls_f';
   g = []; H = [];
   gklsType = 1;
   [x_opt, f_opt, rho, idx_glob] = Tgkls('init',nf,n,nmin, f_star, r_star, rho_star, x_L, x_U);
   Name = sprintf('GKLS-ND %03d',nf);
elseif P >= 101 & P <= 200
   % D problems, f(x) and g(x) are available
   [x_opt, f_opt, rho, idx_glob] = Tgkls('init',nf-100,n,nmin, f_star, r_star, rho_star, x_L, x_U);
   f = 'gkls_f';
   g = 'gkls_g';
   H = [];
   gklsType = 2;
   Name = sprintf('GKLS-D %03d',nf-100);
elseif P >= 201 & P <= 300
   % D2 problems, 2 levels of derivatives available
   [x_opt, f_opt, rho, idx_glob] = Tgkls('init',nf-200,n,nmin, f_star, r_star, rho_star, x_L, x_U);
   f = 'gkls_f';
   g = 'gkls_g';
   H = 'gkls_H';   
   gklsType = 3;
   Name = sprintf('GKLS-D2 %03d',nf-200);
else
   error('gkls_prob: Illegal problem number');
end

% Tgkls returns optimal points as columns in x_opt
x_opt = x_opt';
nLocal = size(x_opt,1);
nGlobal = length(idx_glob);

% Set x_0 to zeros
x_0=zeros(length(x_L),1);

Prob = conAssign(f, g, H, [], x_L,...
    x_U, Name, x_0, [], [], [], [], [], [], [],...
    [], [], [], [], [], [], f_opt, x_opt);
Prob.MIP.nGlobal = nGlobal;
Prob.MIP.nLocal = nLocal;
Prob.gklsType = gklsType;

% MODIFICATION LOG:
%
% 050105 ango Wrote file
% 080603 med  Switched to conAssign, cleaned