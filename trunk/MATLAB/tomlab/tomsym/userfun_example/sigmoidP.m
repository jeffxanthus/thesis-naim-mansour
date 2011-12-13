function P = sigmoidP(f,xp) %#ok
% sigmoidP - Pattern to the logistic sigmoid function.
%
% This function is used by tomSym's "pattern" function, to compute the
% sparsity pattern of the "sigmoid" function.

% This function is provided as an example of how to define your own
% symbolic functions in tomSym.

% In the "pattern" command, it is not allowed to use any information
% symbolic arguments, except their size and their pattern.
% We can call pattern(x) to find the sparsity pattern of x, but since the
% sigmoid function maps zeroes to nonzero, that information is of no value.

P = true(size(xp));

