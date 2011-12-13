% function Prob = modify_c_L(Prob, c_L, idx)
%
% PURPOSE:    Modify lower bounds for nonlinear constraints. If idx is not
%             given c_L will be replaced
%
% INPUTS:
%
% Prob        Existing TOMLAB problem
% c_L         New lower bounds for the nonlinear constraints
% idx         Indices for the modified constraint bounds (optional)
%
% OUTPUTS:
%
% Prob        Modified TOMLAB problem

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Aug 15, 2006.   Last modified Aug 15, 2006.

function Prob = modify_c_L(Prob, c_L, idx)

if nargin < 2
    error('modify_c_L requires at least 2 inputs');
end

if nargin == 3
    if length(idx) ~= length(c_L)
        error('Length of idx does not match length of c_L');
    end
    if max(idx) > Prob.mNonLin
        error('Indices exceed number of nonlinear constraints');
    end
    if min(idx) < 1
        error('Values in idx are lower than 1');
    end
end

if nargin == 2
    if length(c_L) ~= Prob.mNonLin
        error('Length of c_L does not match current length of c_L');
    end
    Prob.c_L = full(double(c_L(:)));
else
    Prob.c_L(idx) = double(c_L(:));
end

if isempty(Prob.c_U)
    Prob.c_U = inf*ones(Prob.mNonLin,1);
end
if any(Prob.c_L>Prob.c_U)
    error('c_L and c_U have crossover values');
end

% MODIFICATION LOG
%
% 060814  med  Written