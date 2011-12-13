% function Prob = add_A(Prob, A, b_L, b_U)
%
% PURPOSE:    Adds linear constraints to an existing problem
%
% INPUTS:
%
% Prob        Existing TOMLAB problem
% A           The additional linear constraints
% b_L         The lower bounds for the new linear constraints
% b_U         The upper bounds for the new linear constraints
%
% OUTPUTS:
%
% Prob        Modified TOMLAB problem

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Aug 15, 2006.   Last modified Aug 15, 2006.

function Prob = add_A(Prob, A, b_L, b_U)

if nargin < 4
    error('add_A requires at least 4 inputs');
end

mLin = max([size(A,1),length(b_L),length(b_U)]);

if size(A,1) ~= mLin | length(b_L) ~= mLin | length(b_U) ~= mLin
    error('Size of A does not match size of bounds');
end

if size(A,2) ~= Prob.N
    error('Incorrect number of columns in A');
end

if any(Prob.b_L>Prob.b_U)
    error('b_L and b_U have crossover values');
end

if issparse(Prob.A) | issparse(A)
    Prob.A = [sparse(Prob.A);sparse(A)];
else
    Prob.A = [Prob.A;A];
end
Prob.b_L = [Prob.b_L;double(b_L(:))];
Prob.b_U = [Prob.b_U;double(b_U(:))];
Prob.mLin = Prob.mLin + mLin;

% MODIFICATION LOG
%
% 060814  med  Written