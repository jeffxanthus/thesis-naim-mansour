% function Prob = modify_b_U(Prob, b_U, idx)
%
% PURPOSE:    Modify upper bounds for linear constraints. If idx is not
%             given b_U will be replaced
%
% INPUTS:
%
% Prob        Existing TOMLAB problem
% b_U         New upper bounds for the linear constraints
% idx         Indices for the modified constraint bounds (optional)
%
% OUTPUTS:
%
% Prob        Modified TOMLAB problem

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Aug 15, 2006.   Last modified Aug 15, 2006.

function Prob = modify_b_U(Prob, b_U, idx)

if nargin < 2
    error('modify_b_U requires at least 2 inputs');
end

if nargin == 3
    if length(idx) ~= length(b_U)
        error('Length of idx does not match length of b_U');
    end
    if max(idx) > Prob.mLin
        error('Indices exceed number of linear constraints');
    end
    if min(idx) < 1
        error('Values in idx are lower than 1');
    end
end

if nargin == 2
    if length(b_U) ~= Prob.mLin
        error('Length of b_U does not match current length of b_U');
    end
    Prob.b_U = full(double(b_U(:)));
else
    Prob.b_U(idx) = double(b_U(:));
end

if isempty(Prob.b_L)
    Prob.b_L = -inf*ones(Prob.mLin,1);
end

if any(Prob.b_L>Prob.b_U)
    error('b_L and b_U have crossover values');
end

% MODIFICATION LOG
%
% 060814  med  Written