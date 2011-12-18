function [sols, iters, activationHist] = OMPDeclip(A, y, N, MclA, theta, maxIters, lambdaStop, solFreq, verbose, OptTol)
% SolveOMP: Orthogonal Matching Pursuit
% Usage
%	[sols, iters, activationHist] = SolveOMP(A, y, N, maxIters, lambdaStop, solFreq, verbose, OptTol,Mp,Mn)
% Input
%	A           Either an explicit nxN matrix, with rank(A) = min(N,n) 
%               by assumption, or a string containing the name of a 
%               function implementing an implicit matrix (see below for 
%               details on the format of the function).
%	y           vector of length n.
%   N           length of solution vector. 
%	maxIters    maximum number of iterations to perform. If not
%               specified, runs to stopping condition (default)
%   lambdaStop  If specified, the algorithm stops when the last coefficient 
%               entered has residual correlation <= lambdaStop. 
%   solFreq     if =0 returns only the final solution, if >0, returns an 
%               array of solutions, one every solFreq iterations (default 0). 
%   verbose     1 to print out detailed progress at each iteration, 0 for
%               no output (default)
%	OptTol      Error tolerance, default 1e-5
% Outputs
%	 sols            solution(s) of OMP
%    iters           number of iterations performed
%    activationHist  Array of indices showing elements entering  
%                    the solution set
% Description
%   SolveOMP is a greedy algorithm to estimate the solution 
%   of the sparse approximation problem
%      min ||x||_0 s.t. A*x = b
%   The implementation implicitly factors the active set matrix A(:,I)
%   using Cholesky updates. 
%   The matrix A can be either an explicit matrix, or an implicit operator
%   implemented as an m-file. If using the implicit form, the user should
%   provide the name of a function of the following format:
%     y = OperatorName(mode, m, n, x, I, dim)
%   This function gets as input a vector x and an index set I, and returns
%   y = A(:,I)*x if mode = 1, or y = A(:,I)'*x if mode = 2. 
%   A is the m by dim implicit matrix implemented by the function. I is a
%   subset of the columns of A, i.e. a subset of 1:dim of length n. x is a
%   vector of length n is mode = 1, or a vector of length m is mode = 2.
global B
global yU


if nargin < 10,
	OptTol = 1e-5;
end
if nargin < 9,
    verbose = 0;
end
if nargin < 8,
    solFreq = 0;
end
if nargin < 7,
    lambdaStop = 0;
end
if nargin < 6,
    maxIters = length(y);
end

explicitA = ~(ischar(A) || isa(A, 'function_handle'));
n = length(y);

% Parameters for linsolve function
% Global variables for linsolve function
global opts opts_tr machPrec
opts.UT = true; 
opts_tr.UT = true; opts_tr.TRANSA = true;
machPrec = 1e-5;

% Initialize
x = zeros(N,1);
k = 1;
R_I = [];
activeSet = [];
sols = [];
res = y;
normy = norm(y);
resnorm = normy;             
done = 0;
offSet=max(abs(y));

while ~done
    if (explicitA)
        corr = A'*res;
    else
        corr = feval(A,2,n,N,res,1:N,N); % = A'*y
    end
    [maxcorr i] = max(abs(corr));
    newIndex = i(1);
    activeSet = [activeSet newIndex];
    %TO CHECK
%     options = optimset('LargeScale','off'); 

%     x = fmins('normax-b.m', x0, options, A, b); ALTERNATIVE
%     MclA(:,activeSet)*(0.1.*(ones(length(activeSet),1)))
%     theta
%     pause
      warning off all;
      if length(activeSet)>maxIters-2%round(maxIters.*0.92)
          options=optimset('Display','final','LargeScale','off','MaxIter',700);
          x0=x(activeSet);
         [x(activeSet),rnorm,res,error]=lsqlin(A(:,activeSet),y,MclA(:,activeSet),theta,[],[],[],[],[x0(1:end-1,1); mean(x0(1:end-1),1)],options);
%          if error==-2
%              theta=0.9.*theta;
%              [x(activeSet),rnorm,res,error]=lsqlin(A(:,activeSet),y,MclA(:,activeSet),theta,[],[],[x(activeSet); mean(x(activeSet))]);
%          end
         if (error==-2 || error==-3 || error==-4 || error==-7)
            x(activeSet)=lsqlin(A(:,activeSet),y);
         end
      else
         x(activeSet)=lsqlin(A(:,activeSet),y);
         disp('Error avoided')
      end
%       J=gradient(y-A(:,activeSet)*x(activeSet));
%       H=J'*J;
%       f=J'*(y-A(:,activeSet)*x(activeSet));
%         %TODO
%       x(activeSet)=quadprog(H,f,MclA(:,activeSet),offSet*ones(length(MclA),1));
%     B=A(:,activeSet);
%     MclB=MclA(:,activeSet);
%     yU=y;
%     options = optimset('Algorithm','interior-point','Display','off');
%     [x(activeSet), fval]=fmincon(@L2Norm,yU(activeSet,:),MclB,offSet*ones(length(MclB),1),[],[],[],[],[],options);

    % Compute new residual
    if (explicitA)
        res = y - A(:,activeSet) * x(activeSet);
    else
        Ax = feval(A,1,n,N,x,1:N,N);
        res = y - Ax;
    end
    resnorm = norm(res);
    
    if ((resnorm <= OptTol*normy) | ((lambdaStop > 0) & (maxcorr <= lambdaStop)))
        length(activeSet)
        done = 1;
    end

    if verbose
        fprintf('Iteration %d: Adding variable %d\n', k, newIndex);
    end

    k = k+1;
    if k >= maxIters
        length(activeSet)
        done = 1;
    end

    if done | ((solFreq > 0) & (~mod(k,solFreq)))
        sols = [sols x];
    end
end

iters = k;
activationHist = activeSet;
clear opts opts_tr machPrec


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, flag] = updateChol(R, n, N, A, explicitA, activeSet, newIndex)
% updateChol: Updates the Cholesky factor R of the matrix 
% A(:,activeSet)'*A(:,activeSet) by adding A(:,newIndex)
% If the candidate column is in the span of the existing 
% active set, R is not updated, and flag is set to 1.

global opts_tr machPrec
flag = 0;

if (explicitA)
    newVec = A(:,newIndex);
else
    e = zeros(N,1);
    e(newIndex) = 1;
    newVec = feval(A,1,n,N,e,1:N,N); 
end

if length(activeSet) == 0,
    R = sqrt(sum(newVec.^2));
else
    if (explicitA)
        p = linsolve(R,A(:,activeSet)'*A(:,newIndex),opts_tr);
    else
        AnewVec = feval(A,2,n,length(activeSet),newVec,activeSet,N);
        p = linsolve(R,AnewVec,opts_tr);
    end
    q = sum(newVec.^2) - sum(p.^2);
    if (q <= machPrec) % Collinear vector
        flag = 1;
    else
        R = [R p; zeros(1, size(R,2)) sqrt(q)];
    end
end

%
% Copyright (c) 2006. Yaakov Tsaig
%  

%
% Part of SparseLab Version:100
% Created Tuesday March 28, 2006
% This is Copyrighted Material
% For Copying permissions see COPYING.m
% Comments? e-mail sparselab@stanford.edu
%
