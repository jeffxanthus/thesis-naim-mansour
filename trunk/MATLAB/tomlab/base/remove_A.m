% function Prob = remove_A(Prob, idx)
%
% PURPOSE:    Removes the linear constraints specified by idx
%
% INPUTS:
%
% Prob        Existing TOMLAB problem
% idx         The row indices to remove in the linear constraints
%
% OUTPUTS:
%
% Prob        Modified TOMLAB problem

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Aug 15, 2006.   Last modified Aug 15, 2006.

function Prob = remove_A(Prob, idx)

if nargin < 2
    error('remove_A requires at least 2 inputs');
end

if max(idx) > Prob.mLin
    error('Indices exceed number of linear constraints');
end

if min(idx) < 1
    error('Values in idx are lower than 1');
end

idx = unique(idx);

Prob.A(idx,:) = [];
Prob.b_L(idx) = [];
Prob.b_U(idx) = [];
Prob.mLin = Prob.mLin-length(idx);

% MODIFICATION LOG
%
% 060814  med  Written