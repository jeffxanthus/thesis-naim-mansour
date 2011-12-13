% function Prob = modify_c(Prob, c, idx)
%
% PURPOSE:    Modify linear objective (LP/QP only)
%
% INPUTS:
%
% Prob        Existing TOMLAB problem
% c           New linear coefficients
% idx         Indices for the modified linear coefficients (optional)
%
% OUTPUTS:
%
% Prob        Modified TOMLAB problem

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Aug 15, 2006.   Last modified Aug 15, 2006.

function Prob = modify_c(Prob, c, idx)

if nargin < 2
    error('modify_c requires at least 2 inputs');
end

if nargin == 3
    if length(idx) ~= length(c)
        error('Length of idx does not match length of c');
    end
    if max(idx) > Prob.N
        error('Indices exceed number of variables');
    end
    if min(idx) < 1
        error('Values in idx are lower than 1');
    end
end

if nargin == 2
    if length(c) ~= Prob.N
        error('Length of c does not match current length of c');
    end
    Prob.QP.c = full(double(c(:)));
else
    Prob.c(idx) = double(c(:));
end

% MODIFICATION LOG
%
% 060814  med  Written