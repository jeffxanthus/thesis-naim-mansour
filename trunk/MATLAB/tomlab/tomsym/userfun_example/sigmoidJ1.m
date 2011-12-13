function J = sigmoidJ1(x)
% sigmoidJ1 - Jacobian to the logistic sigmoid function.
%
% This function is used by tomSym's "derivative" function, to compute the
% Jacobian matrix of the "sigmoid" function.

% This function is provided as an example of how to define your own
% symbolic functions in tomSym.

% The sigmoid function has the derivative: sigmoid(x)*(1-sigmoid(x))

% The Jacobian of an element-wise operation is a diagonal matrix, so we
% use the "setdiag command"

sx = sigmoid(x);

J = setdiag(sx.*(1-sx));


