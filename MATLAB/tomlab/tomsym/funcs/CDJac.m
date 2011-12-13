% CDJac - Numerical approximation of the Jacobian via central differences
%
% J = CDJac(nOut, fun, nIn, h, pattern, colidx, ...) computes the Jacobian
% matrix of FUN with respect to the m:th input argument.
%
% Indata:
% * nOut    - Which output argument to take from FUN (1 = the first one).
% * fun     - A function name or handle, to be called via feval
% * nIn     - Which input argument to compute Jacobian with respect to.
% * h       - The step length to use for central differences (should either
%             be a scalar, or a matrix of the same size as the input.)
% * pattern - The sparsity pattern of the Jacobian
% * colidx  - Indexes, as returned by findpatt(pattern)
% * ...     - All the input arguments to fun...
%
% Example: CDJac(1, 'cos', 1, 1e-7, 1, 1, 0.3) computes the derivative of
% cos(x) at x=0.3. (The true answer is -sin(0.3), but CDJac will only
% return a numerical approximation of the true value of the derivative.)
%
% It is of vital importance that the step length h is properly choosen. An
% incorrect step length may lead to completely incorrect results. If h is
% not given, then h=sqrt(eps)*(1+abs(x)) is used.
%
% See also: FDJac
