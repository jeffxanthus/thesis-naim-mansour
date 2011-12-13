function c = sos2(a)
% tomSym/sos2 - overloaded function
%
% sos2(A) creates a special ordered set condition:
%
% At most 2 elements of a vector A can be strictly positive, the rest being
% zero. Nonzero elements must be consecutive.
%
% For a matrix A, each column forms a special ordered set.
%
% sos2 conditions require the CPLEX or the GUROBI solver.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

if isempty(a)
    c = {};
    return
elseif numel(a)==1
    c = {a>=0};
    return
end

if size(a,1)==1
    a = a'; % Row vector -> column vector
end
    
if size(a,1)==2
    c = {a>=0};
    return
end

c = tomSym(mfilename,size(a,1),size(a,2),a);
