function [f,c] = rewriteV(f,c)
% rewriteV - Rewrite an optimization problem to avoid sharp corners.
%
% [f,c] = rewriteV(f,c) rewrites an optimization problem, characterized by
% an objective function f, and a constraint set c, so as to avoid sharp
% corners in the objective funcition.
%
% f must be a scalar tomSym object.
% c must ba a cell array of tomSym constraints.
%
% Optimization alogrithms typically search for points where the derivative
% of f is zero. If the objective function contains functions like "abs"
% that have sharp corners, then such algorithms may fail to converge, or
% converge very slowly. RewriteV attemtps to rewrite the problem introducing
% extra unknowns and constraints, and eliminating the sharp corners.

% Only the overloaded method has any effect.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-05 by rutquist for TOMLAB release 7.7

[f,c] = rewriteV(tomSym(f),c);
