% function Hx = HxFunc(x, nState, Prob)
%
% Compute H*x for quadratic problem
%
% Input:
% x      Point x where H*x is evaluated 
% nState = 1 if first call
% Prob   Problem structure 
%
% Output:
% Hx     H*x, where H is the matrix in the QP problem (x'  H x) 
%
% H is stored in Prob.QP.F if using the Tomlab format

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Feb 27, 2002.  Last modified Feb 27, 2002.

function Hx = HxFunc(x, nState, Prob)

if isempty(Prob.QP.F)
   % LP Problem
   Hx = zeros(length(x),1);
else
   % QP Problem
   Hx = Prob.QP.F*x;
end