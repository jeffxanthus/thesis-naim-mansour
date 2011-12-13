function y = sos1(a)
% sos1 - Test if a vectors form Special Ordered Sets of Type One 
%
% sos1(a) returns true if at most one member of the vector "a" is strictly 
% positive, and the rest are zero.
%
% If a is a matrix, then sos1 is applied to each column, and a row vector
% of boolean values is returned.
%
% This function mainly exists so that it can be overloaded for tomSym,
% where it is handled efficiently by the CPLEX and GUROBI solvers. When 
% sos1(a) is used as a tomSym constraint, it is mathematically
% equivalent to {a>=0, sum(a>0)<=1}

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

y = all(a>=0) & sum(a>0)<=1;
