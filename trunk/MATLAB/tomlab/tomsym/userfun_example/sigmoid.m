function y = sigmoid(x)
% sigmoid - The logistic sigmoid function.
%
% y = sigmoid(x) computes the logistic function, or "smoothed Heaviside",
% which is the most common sigmoid function. It is uesful in statistics,
% economics, artificial neural networks, and other applications.
%
% It is defined by the equation:
%
%   y = 1./(1+exp(-x)) = 0.5*(1+tanh(0.5*x))

% This function is provided as an example of how to define your own
% symbolic functions in tomSym. This function, would not need to be
% overloded, because it can already accept a tomSym object as input.
% However, by overloading it, it is possible to use another formula for the
% derivative than the one that would result from tomSym's "derivative"
% command.

y = 1./(1+exp(-x));

