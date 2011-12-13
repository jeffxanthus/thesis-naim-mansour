% function Prob = replace_A(Prob, A, b_L, b_U)
%
% PURPOSE:    Replaces the linear constraints
%
% INPUTS:
%
% Prob        Existing TOMLAB problem
% A           New linear constraints
% b_L         Lower bounds for linear constraints
% b_U         Upper bounds for linear constraints
%
% OUTPUTS:
%
% Prob        Modified TOMLAB problem

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Aug 15, 2006.   Last modified Aug 15, 2006.

function Prob = replace_A(Prob, A, b_L, b_U)

if nargin < 3
    error('replace_A requires at least 3 inputs');
end

if nargin == 3
    b_U = [];
end

[m,mN]=size(A);
if mN ~= Prob.N
    fprintf('Number of variables %d\n',Prob.N);
    fprintf('Number of columns in linear constraint matrix A %d\n',mN);
    fprintf('These lengths should be the same, check input.\n');
    error('Illegal number of columns in A');
end
Prob.mLin = m;
Prob.A = A;
if isempty(b_L)
    Prob.b_L=-Inf*ones(m,1);
elseif m==length(b_L)
    Prob.b_L=full(double(b_L(:)));
else
    fprintf('Length of b_L %d, Rows in A are %d\n',length(b_L),m);
    error('Illegal length of b_L');
end
if isempty(b_U)
    Prob.b_U=Inf*ones(m,1);
elseif m==length(b_U)
    Prob.b_U=full(double(b_U(:)));
else
    fprintf('Length of b_U %d, Rows in A are %d\n',length(b_U),m);
    error('Illegal length of b_U');
end
if any(Prob.b_L>Prob.b_U)
    error('b_L and b_U have crossover values');
end

% MODIFICATION LOG
%
% 060814  med  Written