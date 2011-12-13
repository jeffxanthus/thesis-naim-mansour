% checkAssign is a verification routine for input data to assign routines
% included with TOMLAB.
%
% Prob = checkAssign(Prob, ...)
%
% -----------------------------------------------------
%
% Syntax of checkAssign:
%
% function Prob = checkAssign(Prob, x_0, x_L, x_U, b_L, b_U, A);
%
% INPUT (One parameter c must always be given)
%
% Prob:        The Problem structure
% x_0:         Starting point x (may be empty)
% x_L:         Lower bounds on x
% x_U:         Upper bounds on x
% b_L:         The lower bounds for the linear constraints
% b_U:         The upper bounds for the linear constraints
% A:           The linear constraint matrix
%
%              b_L, b_U, x_L, x_U must either be empty or of full length

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2006-2007 by Tomlab Optimization Inc., Sweden. $Release: 5.7.0$
% Written Dec 13, 2006. Last modified Dec 13, 2006.

function Prob = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A)

if isempty(x_0)
    Prob.x_0 = zeros(n,1);
elseif length(x_0) ~= n
    fprintf('Length of x_0 %d, should be %d\n',length(x_0),n);
    error('Illegal length of x_0');
else
    Prob.x_0 = full(double(x_0(:)));
end

if isempty(x_L)
    Prob.x_L = -Inf*ones(n,1);
elseif length(x_L) ~= n
    fprintf('Length of x_L %d, should be %d\n',length(x_L),n);
    error('Illegal length of x_L');
else
    Prob.x_L = full(double(x_L(:)));
end

if isempty(x_U)
    Prob.x_U = Inf*ones(n,1);
elseif length(x_U) ~= n
    fprintf('Length of x_U %d, should be %d\n',length(x_U),n);
    error('Illegal length of x_U');
else
    Prob.x_U = full(double(x_U(:)));
end
if any(Prob.x_L>Prob.x_U)
    error('x_L and x_U have crossover values');
end

if ~isempty(A)
    Prob.A=A;
    [m,mN]=size(A);
    if mN ~= n
        fprintf('Number of variables %d\n',n);
        fprintf('Number of columns in linear constraint matrix A %d\n',mN);
        fprintf('These lengths should be the same, check input!!!\n');
        error('Illegal number of columns in A')
    end
    Prob.mLin = m;
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
else
    Prob.mLin = 0;
end

% MODIFICATION LOG
%
% 061213  med  Written