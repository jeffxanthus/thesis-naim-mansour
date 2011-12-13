% function dc = oqnlp_dc(x, Prob, varargin)
%
% TOMLAB interface routine to avoid numerical Jacobian differentiation of
% integer variables in MINLP problems

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Nov 13, 2004.   Last modified Nov 13, 2004.

function dc = oqnlp_dc(x, Prob, varargin)

IntVars    = Prob.MIP.IntVars;
m          = Prob.mNonLin;
if m > 0
    dc            = NaN*Prob.ConsPattern;
    dc(:,IntVars) = 0;
else
    dc            = [];
end

% MODIFICATION LOG
%
% 041113  med  Written