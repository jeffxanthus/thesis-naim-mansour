% CPLEX MEX-interface internal callback routine
%
% Makes the cpxCBInfo vector global

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 11.0.0$
% Written Aug. 8, 2002 Last modified Mar 22, 2005

function cpx2cbinfo(a)

global cpxCBInfo

cpxCBInfo = a;
