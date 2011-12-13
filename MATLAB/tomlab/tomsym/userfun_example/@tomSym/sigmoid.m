function y = sigmoid(x)
% tomSym/sigmoid - Overloaded logistic sigmoid function.

% This function is provided as an example of how to define your own
% symbolic functions in tomSym. The overloaded function only needs to
% return a symbolic object of the correct size.

y = tomSym('sigmoid',size(x,1),size(x,2),x);


