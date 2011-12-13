% effiFrontEx.m:
%
% Example for QP optimization with nonsmooth constraints:
%
% bl <= sum(abs(x-x_s) <= bu

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 18, 2004.   Last modified Sep 18, 2005.

clear all
load portfolioAbsConExMat % Load test data

% Assign problem without additional constraints
Prob = qpAssign(sparse(H+H'),f,sparse(A),[],b,lb,ub);

% Add abs constraints.
if 1
   Prob = portfolioAbsCon(Prob, [], 0.35, xi);
else
   % It works for 0.07, but not for 0.06
   Prob = portfolioAbsCon(Prob, [], 0.07, xi);
end

% Choose a sparse QP solver depending on license as the problem is sparse
Solver = GetSolver('qp', 1, 1);

% Solve the problem
Result = tomRun(Solver,Prob,1);
%Result = tomRun('snopt',Prob,1);
%Result = tomRun('cplex',Prob,1);

format long

% Check the solution
N = Prob.N/2;
x_k = Result.x_k;

% Round to 0 any variables close to 0, use limit 1E-10
ix = find(x_k < 1E-10);
x_k(ix) = 0;

z=double(x_k(1:N) > 0)+ double(x_k(N+1:2*N) > 0);
disp('check for both positive and negative values');
idx = find(z > 1);
if ~isempty(idx)
    idx
end

z = x_k(1:N) - x_k(N+1:2*N) + xi(:);
disp('check for invalid linear constraints');
ix=find( A*z > b+1E-5);
if ~isempty(ix)
   resid = A*z-b;
   resid(ix)
end

disp('check for invalid bounds');
ub = ub(:);
lb = lb(:);

iub = find( z > ub(:)+1E-5);
if ~isempty(iub)
   [iub lb(iub) z(iub) ub(iub)]
end
ilb = find( z < lb(:)-1E-5);
if ~isempty(iub)
   [ilb lb(ilb) z(ilb) ub(ilb)]
end

% MODIFICATION LOG:
%
% 031014 hkh  Written
% 041221 hkh  Revised. Round down small x_k to 0, such x_k if barrier QP
% 050919 med  Constant tranferred to Prob.fConstant
