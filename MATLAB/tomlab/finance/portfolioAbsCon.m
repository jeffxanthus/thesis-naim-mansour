% function [Prob, fCons] = portfolioAbsCon(Prob, bl, bu, x_s)
%
% For a given quadratic program (QP), adds nonsmooth constraint
%
% bl <= sum(abs(x-x_s)) <= bu
%
% and reformulates the problem to a smooth problem with twice as many variables
%
% Output:
%
% Prob            Tomlab problem structure with reformulated problem
%
% Prob.fConstant  Constant factor in objective function after transformation

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2006 by Tomlab Optimization Inc., $Release: 5.4.0$
% Written Oct 14, 2003. Last modified Aug 8, 2006.

function Prob = portfolioAbsCon(Prob, bl, bu, x_s)

N = Prob.N;

if nargin < 4
   x_s = [];
   if nargin < 3
      error('PortfolioAbsCon needs at least three parameters');
   end
end

if isempty(x_s), x_s = zeros(N,1); end
if isempty(bl),  bl  = -inf; end
if isempty(bu),  bu  =  inf; end
x_s = x_s(:);

m = size(Prob.A,1);
b_L = Prob.b_L(:);
b_U = Prob.b_U(:);
if isempty(b_L), b_L = -inf*ones(m,1); end
if isempty(b_U), b_U =  inf*ones(m,1); end

x_L = Prob.x_L(:);
x_U = Prob.x_U(:);
if isempty(x_L), x_L = -inf*ones(N,1); end
if isempty(x_U), x_U =  inf*ones(N,1); end

c = DefPar(Prob.QP,'c',[]);

if isempty(c)
   Prob.fConstant = 0;
   c     = zeros(N,1);
else
   Prob.fConstant = c(:)'*x_s;
end

if ~isempty(Prob.QP.F)
   c = c(:)+2*Prob.QP.F*x_s;
   Prob.fConstant = Prob.fConstant + x_s'*Prob.QP.F*x_s;
   Prob.QP.F = [Prob.QP.F -Prob.QP.F ; -Prob.QP.F Prob.QP.F];
end

if ~isempty(c)
   Prob.QP.c = [c(:) ; -c(:)];
end

if ~isempty(Prob.A)
   Ax = Prob.A*x_s;
   Prob.A = [sparse(Prob.A) -sparse(Prob.A);sparse(ones(1,2*N));speye(N),-speye(N)];
   Prob.b_L = [b_L-Ax ; bl; x_L - x_s ];
   Prob.b_U = [b_U-Ax ; bu; x_U - x_s ];
else
   Prob.A   = [sparse(ones(1,2*N));speye(N),-speye(N)];
   Prob.b_L = [bl; x_L + x_s ];
   Prob.b_U = [bu; x_U + x_s ];
end

Prob.x_U = inf*ones(2*N,1);
Prob.x_L = zeros(2*N,1);

if ~isempty(Prob.x_0)
   Prob.x_0 = [max(0,Prob.x_0(:)-x_s);-min(0,Prob.x_0(:)-x_s)];
end

Prob.mLin = max(length(Prob.b_L), length(Prob.b_U));
Prob.N = 2*N;

% MODIFICATION LOG:
%
% 031014 hkh  Written
% 041221 hhk  Mapping of x_0 fixed
% 050117 med  mlint revision
% 050919 med  Make sure all matrices are sparse directly
% 050919 med  fCons transferred to Prob.fConstant
% 060808 med  Updated help