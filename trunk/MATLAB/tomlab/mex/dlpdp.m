%    dlpdp solver
% 
%    dlpdp solves the least projected distance problem. 
%     
%    Given dense matrices G and H of respective
%    dimensions m by N1 and m by N2, and vectors y of
%    length m.  This subroutine solves the
%    linearly constrained least squares problem formulated as:
%
%             min   || w || subject to
%             w,z
%
%                 G w + H z >= y    Inequality constraints
%
%    Any number of rows in G and H are permitted.
%
%
% function [x, wNorm, mode] = dlpdp ( ...
%     G, H, y, Scale, ColScale, BlowUp, RankTol)
%
% Input: (at least 2 input parameters needed)
%   G         m  x n1 dense matrix
%   H         m  x n2 dense matrix
%   y         m  x 1 dense vector
%   Scale     If > 0 scale all nonzero values in G and H. Default 1
%   ColScale  If nonempty, n x 1 dense vector with diagonal scaling of columns
%   RankTol   Rank determination tolerance
%   BlowUp    Blow-up parameter
%
% Output:
%   x         Solution x
%   wNorm     ||w||
%   mode      Exit status:      
%             0  The solution was successfully obtained.
%             1  The inequalities are inconsistent.
%

%   Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
%   Copyright (c) 2000-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
%   Written Jul 17, 2009. Last modified Jul 17, 2009.
%# mex

function [x, wNorm, mode] = dlpdp ( ...
    G, H, y, Scale, ColScale, RankTol, BlowUp)

help dlpdp;