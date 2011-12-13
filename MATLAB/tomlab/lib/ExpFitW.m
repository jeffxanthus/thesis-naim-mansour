% function W = ExpFitW(x, r, Prob)
%
% Weighting routine used in exponential fitting problems
% 
% The name of this routine is defined in Prob.LS.weightY
%
% The routine is called from the gateway routine nlp_r
%
% INPUT:
% x       Current iterate x, where the residual is vector function r=r(x).
% r       Residual
% Prob    TOMLAB problem structure
%
% OUTPUT:
% W       Weight vector or matrix

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2011 by Tomlab Optimization Inc., Sweden. $Release: 7.8.0$
% Written Oct 24, 1998.  Last modified July 22, 2011.

function W = ExpFitW(x, r, Prob)

[SepAlg, p, wType] = expGet1(Prob);

y     = Prob.LS.y;
m     = size(y,1);
W     = ones(m,1);

if wType==1
   ix    = find(y~=0 & ~isinf(y));
   W(ix) = 1./abs(y(ix));
elseif wType==2  
   r     = r(:);
   yMod  = r+y;
   ix    = find(yMod~=0);
   W(ix) = 1./abs(yMod(ix));
elseif wType==3  
   r     = r(:);
   yMod  = r+y;
   z     = max(yMod,y);
   ix    = find(z~=0);
   W(ix) = 1./abs(z(ix));
end

% MODIFICATION LOG:
%
% 981025  hkh  r not accessed for wType==0,1. Avoid using. 
% 981108  hkh  Remove wrong comments 
% 981129  hkh  YMod should be yMod on two places.
% 981207  hkh  Must make fix to avoid weights to be 0 when plotUtil calls.
% 110723  hkh  Revised. Avoid division with 0 for all wType
