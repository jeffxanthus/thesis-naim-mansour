% expSet2 Set values of four parameters for exponential fitting in Prob.ExpFit
%
% The parameters are used to find starting values
%
% function Prob = expSet2(Prob, dType,geoType,qType,sigType); 
%
% INPUT:
%  Prob     Problem structure
%  dType
%  geoType
%  qType
%  sigType
% OUTPUT:
%  Prob     Problem structure, changed entries in Prob.ExpFit
           
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Feb 19, 1998.  Last modified June 21, 1999.

function Prob = expSet2(Prob, dType,geoType,qType,sigType)

if nargin > 1
   Prob.ExpFit.dType=dType;
   if nargin > 2
      Prob.ExpFit.geoType=geoType;
      if nargin > 3
         Prob.ExpFit.qType=qType;
         if nargin > 4
            Prob.ExpFit.sigType=sigType;
         end
      end
   end
end