% function [f,c] = glcd_fc(x,fname,cname,varargin)
%
% TOMLAB gateway routine for the computation of objective function and
% nonlinear constraint values for the glcDirect solver.

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2007-2007 by Tomlab Optimization Inc., $Release: 6.0.0$
% Written Sep 5, 2007.   Last modified Sep 6, 2007.

function [f,c]=glcd_fc(x,fname,cname,varargin)
f=feval(fname,x,varargin{:});
c=feval(cname,x,varargin{:});

% MODIFICATION LOG
%
% 070905 ango Wrote file