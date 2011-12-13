function y = tomNorm(v,p)
%tomNorm - A vector norm function that is efficient for optimization
%
% tomNorm(V,P) >= sum(abs(V).^P)^(1/P)
%
% tomNorm(V,2) or tomNorm(V) gives the 2-norm of V
% tomNorm(V,1) gives the 1-norm of V.
% tomNorm(V,inf) is the infinity norm of V.
%
% V must be a real vector (not a matrix). P must be greater than zero.
%
% tomNorm is similar to Matlab's built-in norm function for vectors,
% but it can only be used in an optimization context.
%
% When tomNorm is used in either the objective of constraints of an
% optimization problem, it introduces extra decision variables and
% constraints so that tomNorm is greater than or equal to the real norm.
%
% If tomNorm is an objective to be minimized (which is typically the case)
% then the solution will be the same as if norm had been used, but
% convergence will be much faster. Similarly, constraints requiring a norm
% to be smaller than a given value will benefit from using tomNorm.
%
% If the objective involves maximizing a norm, then tomNorm cannot be used
% tp compute that norm. Similarly, if a constraint requires a norm to be 
% greater than a certain value, then tomNorm cannot be used.
%
% The reason why tomNorm gives better performance is that it avoids using
% the absolute value (abs) function, which is nonlinear with a
% discontinuous derivative. In many cases, tomNorm yields a linear
% programming where norm would yield a nonlinear programming with poor
% convergence.
%
% For L2 norms, where the abs function is not needed,
% tomNorm is the same as norm.


if nargin<2
    p = 2;
end

if numel(v)~=length(v)
    error('V must be a vector');
end

if p<=0
    error('P must be greater than zero.')
elseif mod(p,2) == 0
    y = sum(v.^p)^(1/p);
elseif isinf(p)
    maxabsv = tom([],1,1);
    y = subjectTo(maxabsv,maxabsv==0,maxabsv>=v,maxabsv>=-v);
else
    absv = tom([],size(v,1),size(v,2));
    y = subjectTo(sum(absv.^p)^(1/p),absv==0,absv>=v,absv>=-v);
end
