% tomSym/derivative - The symbolic derivative of a tomSym object.
%
% df_dx    = derivative(f,x) computes the derivative of a tomSym object f
% with respect to a symbol x.
%
% d2f_dxdy = derivative(f,x,y,..) computes the second (third, ...)
% derivative of f with respect to the symbols x, y, ...
%
% Both f and x can be matrices.
%
% The symbol x can either be a tomSym symbol, or a concatenation of tomSym
% symbols.
%
% There exist several conventions for how the elements are arranged in the
% derivative of a matrix. This function uses the convention that vec(df) =
% df_dx*vec(dx). This means that if size(f) = [m n] and size(x) = [p q]
% then size(df_dx) = [m*n, p*q]. Thus, the derivative of a matrix or vector
% with respect to a scalar is a column vector, the derivative of a scalar
% is a row vector, and the derivative of any matrix with respect to itself
% is an identity matrix.
%
% Examples:
% - If f and x are vectors, then J = derivative(f,x) computes the Jacobian
%   matrix.
% - If f is scalar and x is a vector, then H = derivative(f,x,x) computes
%   the Hessian matrix.
%
% For functions that tomSym/derivative is unfamiliar with, it assumes that
% there also exists a derivative function for each argument. For example,
% when computing the derivative of a function "userfun", it is assumed that
% "userfunJ1" gives the derivative of userfun with respect to its first
% input argument.
%
% Reference: Brookes, M., "The Matrix Reference Manual"
%            http://www.ee.ic.ac.uk/hp/staff/dmb/matrix/intro.html
%
% See also: derivatives
