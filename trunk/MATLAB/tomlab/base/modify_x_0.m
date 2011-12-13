% function Prob = modify_x_0(Prob, x_0, idx)
%
% PURPOSE:    Modify starting point. If x_0 is outside the bounds an error
%             will be returned. If idx is not given x_0 will be replaced.
%
% INPUTS:
%
% Prob        Existing TOMLAB problem
% x_0         New starting points
% idx         Indices for the modified starting points (optional)
%
% OUTPUTS:
%
% Prob        Modified TOMLAB problem

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Aug 15, 2006.   Last modified Aug 15, 2006.

function Prob = modify_x_0(Prob, x_0, idx)

if nargin < 2
    error('modify_x_0 requires at least 2 inputs');
end

if nargin == 3
    if length(idx) ~= length(x_0)
        error('Length of idx does not match length of x_0');
    end
    if max(idx) > Prob.N
        error('Indices exceed number of variables');
    end
    if min(idx) < 1
        error('Values in idx are lower than 1');
    end
end

if nargin == 2
    if length(x_0) ~= Prob.N
        error('Length of x_0 does not match current length of x_0');
    end
    Prob.x_0 = full(double(x_0(:)));
else
    Prob.x_0(idx) = double(x_0(:));
end

if isempty(Prob.x_L)
    Prob.x_L = -inf*ones(Prob.N,1);
end

if isempty(Prob.x_U)
    Prob.x_U = inf*ones(Prob.N,1);
end

if any(Prob.x_0>Prob.x_U)
    error('x_0 has values greater than x_U');
end

if any(Prob.x_0<Prob.x_L)
    error('x_0 has values lower than x_L');
end

% MODIFICATION LOG
%
% 060814  med  Written