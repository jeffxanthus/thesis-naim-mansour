% function c = glcd_c(x, func, varargin)
%
% TOMLAB gateway routine for the computation of nonlinear constraint values
% c(x) for the glcDirect solver.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2007-2007 by Tomlab Optimization Inc., $Release: 6.0.0$
% Written Sep 5, 2007.   Last modified Sep 6, 2007.

function c=glcd_c(x,func,varargin)
c=feval(func,x,varargin{:});

% MODIFICATION LOG
%
% 070905 ango Wrote file