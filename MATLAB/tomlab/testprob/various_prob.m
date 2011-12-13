% various_prob:
%
% Various problems
%
% function [probList, Prob] = various_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function [probList, Prob] = various_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'L1LinSolve test example' ...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return;
end

if P == 1
    Name='L1LinSolve test example';
    n = 6;
    x_L = -10*ones(n,1);
    x_U =  10*ones(n,1);
    x_0 = (x_L + x_U) / 2;
    C = spdiags([1 2 3 4 5 6]', 0, n, n);
    y = 1.5*ones(n,1);
    A = [1 1 0 0 0 0]; b_L = 1; b_U = 1;
    Prob = llsAssign(C, y, x_L, x_U, Name, x_0, ...
        [], [], [], ...
        A, b_L, b_U);
    Prob.LS.damp = 1;
    Prob.LS.L = spdiags(ones(6,1)*0.01, 0, 6, 6);
    Prob.SolverL1 = 'lpSimplex';
else
    error('various_prob: Illegal problem number')
end
