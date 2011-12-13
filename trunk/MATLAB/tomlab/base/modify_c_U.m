% function Prob = modify_c_U(Prob, c_U, idx)
%
% PURPOSE:    Modify upper bounds for nonlinear constraints. If idx is not
%             given c_U will be replaced
%
% INPUTS:
%
% Prob        Existing TOMLAB problem
% c_U         New upper bounds for the nonlinear constraints
% idx         Indices for the modified constraint bounds (optional)
%
% OUTPUTS:
%
% Prob        Modified TOMLAB problem

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Aug 15, 2006.   Last modified Aug 15, 2006.

function Prob = modify_c_U(Prob, c_U, idx)

if nargin < 2
    error('modify_c_U requires at least 2 inputs');
end

if nargin == 3
    if length(idx) ~= length(c_U)
        error('Length of idx does not match length of c_U');
    end
    if max(idx) > Prob.mNonLin
        error('Indices exceed number of nonlinear constraints');
    end
    if min(idx) < 1
        error('Values in idx are lower than 1');
    end
end

if nargin == 2
    if length(c_U) ~= Prob.mNonLin
        error('Length of c_U does not match current length of c_U');
    end
    Prob.c_U = full(double(c_U(:)));
else
    Prob.c_U(idx) = double(c_U(:));
end

if isempty(Prob.c_L)
    Prob.c_L = -inf*ones(Prob.mNonLin,1);
end
if any(Prob.c_L>Prob.c_U)
    error('c_L and c_U have crossover values');
end

% MODIFICATION LOG
%
% 060814  med  Written