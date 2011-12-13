% FDJacR - Numerical approximation of the Jacobian
%
% J = FDJac(FUN, info, m, varargin) computes the Jacobian matrix of FUN
% with respect to the m:th input argument.
%
% Implementation based on the algorithm FD in Practical Optimization, page 343.
%
% Based on the FDJac function, but modified to allow reentrant calls and
% work with tomSym.
% FDJacR stores data between evaluations in a persistent variable linked to
% the id field of the info structure. This field is set automatically by
% the WRAP functio.)
%
% See also: wrap
