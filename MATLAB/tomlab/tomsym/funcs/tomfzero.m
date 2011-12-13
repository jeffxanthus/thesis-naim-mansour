function x = tomfzero(fun,x0)
% tomfzero - tomSym varation of fsolve/fzero
%
% Don't use fsolve/fzero is used within the objective function of an 
% optimization problem - Use tomfzero instead.
%
% X = TOMFZERO(FUN,X0) finds X so that FUN(X) equals zero, using X0 as a
% starting guess. X0 can be a scalar, vector or matrix. The returned X
% will have the same size as X0.
%
% FUN should be a function handle, a function name, or a tomSym expression.
%
% If FUN is a tomSym expression, then X0 must be a tomSym equality
% wose left-hand side is the variable to be solved for.
%
% TOMFZERO can combine the equation with other simultaneous constraints and 
% an optimization objective. This can make the solution process orders of 
% magnitude faster compared to using fzero recursively.
%
% Example: Both of the below lines compute f so that sin(f) == y.
%   f = tomfzero(@(x) sin(x)-y, 0)
%   toms x; f = tomfzero(sin(x)==y, x==0)
%
% See also: subjectTo

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

if isa(fun,'tomSym')
    if ~tomCmp(x0,'eq')
        error('Initial guess must be a tomSym equation.');
    end
    xsym = operand(1,x0);
    if ~tomCmp(xsym,'tom');
        error('tomfzero only solves for one unknown');
    end
    if ~tomCmp(fun,'eq')
        fun = (fun==0);
    end
else
    xsym = tom([],size(x0,1),size(x0,2));
    fun = (feval(fun,xsym) == 0);
    x0 = (xsym == x0);
end

x = subjectTo(xsym,x0,fun);


