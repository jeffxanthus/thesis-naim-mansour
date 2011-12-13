% cpx2mat reads (X)MPS, LP and SAV files. The file is converted to matrices
% and vectors made available in MATLAB.
%
% USAGE:
%
%  [F, c, A, b_L, b_U, x_L, x_U, IntVars] = cpx2mat(Name,PriLev);
%
% INPUT:
%
%  Name          Name of the MPS, LP or SAV file with extension.
%
%  PriLev        Print level of cpx2mat. Set to 0 to have it silent, 1 to
%                print warnings, and 2 to print debug information.
%
% OUTPUT:
%
%  F             The quadratic term matrix. Empty for non-QP problems.
%  c             The linear term vector.
%
%  A             The constraint matrix.
%  b_L           The lower bounds of the constraints.
%  b_U           The upper bounds of the constraints.
%
%  x_L           The lower box bounds of x.
%  x_U           The upper box bounds of x.
%
%  IntVars       Logical vector describing what variables that are integer
%                or binary variables. Empty if the problem is not a mixed
%                integer problem.

% Fredrik Hellman, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2006 by Tomlab Optimization Inc., $Release: 11.2.0$
% Written Nov 10, 2004. Last modified Apr 30, 2008.

function [F,c,A,b_L,b_U,x_L,x_U,IntVars] = cpx2mat(Name,PriLev)

if nargin < 2
    PriLev=0;
    if nargin < 1
        error('cpx2mat needs an input file.');
    end
end

[F,c,A,b_L,b_U,x_L,x_U,IntVars] = cpx2matmex(Name,PriLev);

% MODIFICATION LOG
%
% 041110 frhe Written. Used xpr2mat.m as template
% 090430 med  LP and SAV files added