function s = getSolution(result, n)
% getSolution - Extract variables from a result retuned by tomRun.
%
% s = getSolution(result) returns a struct where each field corresponds
% to one of the symbols used in the tomSym problem formulation.
%
% s = getSolution(result, n) returns the n:th solution when the solver
% returned multiple solutions (as is the case with multiMin).
% The index n must be and integer between 1 and size(result.x_k,2)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-09-16 by rutquist for TOMLAB release 7.7

if ~isfield(result.Prob,'tomSym')
    error('The result does not seem to originate from a tomSym problem.');
end

idx = result.Prob.tomSym.idx;
fn = fieldnames(idx);
s = struct;

if nargin < 2
    n = 1;
end

if isfield(result.Prob.tomSym,'xScale')
    xs = result.Prob.tomSym.xScale.*result.x_k(:,n);
else
    xs = result.x_k(:,n);
end

for i=1:length(fn)
    s.(fn{i}) = reshape(xs(idx.(fn{i})),size(idx.(fn{i})));
end
