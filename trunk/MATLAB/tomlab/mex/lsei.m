%    lsei solver
% 
%    lsei solves a linearly constrained least squares
%    problem with both equality and inequality constraints. 
%     
%    May return a covariance matrix of the solution parameters.
%
%    Given dense matrices E, A and G of respective
%    dimensions mE by N, mA by N and mG by N, and vectors f, b and h of
%    respective lengths ME, MA and MG.  This subroutine solves the
%    linearly constrained least squares problem formulated as:
%
%             min   || A x  - b || subject to
%
%                      E x  = f    Equations to be exactly satisfied
%                      G x >= h    Inequality constraints
%
%    In case the equality constraints cannot be satisfied, a
%    generalized inverse solution residual vector length is obtained
%    for f-Ex.  This is the minimal length possible for f-Ex.
%
%    Any number of rows in A, E and G are permitted.
%
%
% function [x, rNormE, rNormL, mode, Cov, eqRank, rlsRank, Iter] = lsei ( ...
%     A, b, E, f, G, h, CompCov, ScaleCov, ColScale, D, epsEQ, epsRank)
%
% Input: (at least 2 input parameters needed)
%   A         mA x n dense matrix
%   b         mA x 1 dense vector
%   E         mE x n dense matrix
%   f         mE x 1 dense vector
%   G         mG x n dense matrix
%   h         mG x 1 dense vector
%   CompCov   If > 0 compute covariance matrix. Default 0
%   ScaleCov  If > 0 scale covariance matrix. Default 1
%   ColScale  If > 0 column scaling of full matrix W. Default 0
%   D         If nonempty, n x 1 dense vector with diagonal scaling of columns
%   epsEQ     Linear equality feasibility tolerance
%   epsRank   Rank tolerance in least squares part
%     
% Output:
%   x         Solution x
%   rNormE    ||Ex-f||
%   rNormL    ||Ax-b||
%   mode      Exit status:      
%             0  Both equality and inequality constraints are compatible 
%                and have been satisfied.
%             1  Equality constraints are contradictory.
%                A generalized inverse solution of Ex=f was used
%                to minimize the residual vector length f-Ex.
%             2  Inequality constraints are contradictory.
%             3  Both equality and inequality constraints are contradictory.
%             4  Usage error occurred. This should not occur.
%
%   Cov       Covariance matrix, if CompCov > 0
%   eqRank    Rank of E matrix part
%   rlsRank   Rank of A matrix part
%   Iter      Number of iterations

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
%   Written Oct 30, 2000. Last modified Sep 10, 2009.
%# mex

function [x, rNormE, rNormL, mode, Cov, eqRank, rlsRank, Iter] = lsei ( ...
    A, b, E, f, G, h, CompCov, ScaleCov, ColScale, D, epsEQ, epsRank)

help lsei;
