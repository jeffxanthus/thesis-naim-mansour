function y = sos2(a)
% sos2 - Test if a vectors form Special Ordered Sets of Type Two
%
% sos2(a) returns true if at most two members of the vector "a" are 
% strictly positive, the rest are zero, and the nonzero elements are
% consecuitive in the vector.
%
% If a is a matrix, then sos2 is applied to each column, and a row vector
% of boolean values is returned.
%
% This function mainly exists so that it can be overloaded for tomSym,
% where it is handled efficiently by the CPLEX and GUROBI solvers. When
% sos2 is used on a tomSym variable it returns a constraint requiring that
% variable to form a special ordered set of type two.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

if size(a,1)==1
    % Row vector
    y = all(a>=0) & sum(a>0)<=2 & sum(diff([0 (a>0)])>0) <= 1;
else
    % Column or matrix
    y = all(a>=0) & sum(a>0)<=2 & sum(diff([zeros(1,size(a,2)); (a>0)])>0) <= 1;    
end
    
