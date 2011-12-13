% expGet2 picks up four parameters for exponential fitting from Prob structure
%
% The parameters are used in finding starting values
%
% function [dType,geoType,qType,sigType] = expGet2(Prob);
%
% INPUT:
%  Prob     Problem structure
%
% OUTPUT:
%  dType
%  geoType
%  qType
%  sigType

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Feb 17, 1998.  Last modified June 21, 1999.

function [dType,geoType,qType,sigType] = expGet2(Prob)

if isempty(Prob)
    dType=  0;
    geoType=0;
    qType=  0;
    sigType=0;
else
    dType=  Prob.ExpFit.dType;
    geoType=Prob.ExpFit.geoType;
    qType=  Prob.ExpFit.qType;
    sigType=Prob.ExpFit.sigType;
end