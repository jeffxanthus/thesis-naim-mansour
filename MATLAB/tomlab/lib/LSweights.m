% LSweights
%
% Example of weighting routine to be used in weighted nonlinear least squares
%
% The name of this routine is set in Prob.LS.weightY
%
% The routine is called from the gateway routine nlp_r
%
% function W = LSweights(x, r, Prob)
%
% INPUT:
% x       Current iterate x, where the residual is vector function r=r(x).
% r       Residual vector, m x 1
% Prob    TOMLAB problem structure
%
% OUTPUT:
% W       Weight vector or matrix

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Oct 18, 1998.  Last modified July 22, 2011.

function W = LSweights(x, r, Prob)

r     = r(:);
m     = length(r);

% Now simply just weight with data
y     = Prob.LS.y;  
W     = ones(m,1);
ix    = find(y~=0 & ~isinf(y));

% Multiplicative weights, W.*r will be applied later
W(ix) = 1./abs(y(ix));

% MODIFICATION LOG
%
% 981207  hkh  Must make fix to avoid weights to be 0 when plotUtil calls.
% 090813  med  mlint check
% 110722  hkh  Revised. W should be multiplicative
