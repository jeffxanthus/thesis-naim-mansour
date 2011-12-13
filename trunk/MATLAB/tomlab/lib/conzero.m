% conzero.m:
%
% Used to find x and y so constraint index equals zero. 
%
% Called from PlotUtil.m
%
% function c = conzero(x, Prob, varargin)
 
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1997-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Nov 16, 1997.   Last modified Sep 30, 2000.

function c = conzero(x, Prob, varargin)

z  = Prob.f0.z;
ix = Prob.f0.ix;

if Prob.f0.var==1
   z(ix(1)) = x;
   z(ix(2)) = Prob.f0.y;
else
   z(ix(1)) = Prob.f0.y;
   z(ix(2)) = x;
end

cc = nlp_cF(z,Prob, varargin{:});
c  = cc(Prob.f0.index);

% MODIFICATION LOG:
%
% 990626  hkh  Avoid feval
% 990909  hkh  Add varargin
% 000306  hkh  Direct call to nlp_cF
% 000930  hkh  Revised for use with Tfzero instead of fzero