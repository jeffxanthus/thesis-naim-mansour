% expSet1 Set values of four parameters for exponential fitting in Prob.ExpFit
%
% function Prob = expSet1(Prob, SepAlg, p, wType, eType, infCR);
%
% INPUT:
%  Prob     Problem structure
%           See userpar.m for a description
%  SepAlg   Flag if Separable NLLS solution
%  p        No of exponential terms
%  wType    Weighting type
%  eType    Type of exponential terms
%  infCR    Information criteria for selection of best number of terms
%
% OUTPUT:
%  Prob     Problem structure, changed entries in Prob.ExpFit
               
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Feb 18, 1998.  Last modified Nov 5, 2000.

function Prob = expSet1(Prob, SepAlg, p, wType, eType, infCR)

if nargin > 1
   Prob.LS.SepAlg=SepAlg;
   if nargin > 2
      Prob.ExpFit.p=p;
      if nargin > 3
         Prob.ExpFit.wType=wType;
         if nargin > 4
            Prob.ExpFit.eType=eType;
            if nargin > 5
               Prob.ExpFit.infCR=infCR;
            end
         end
      end
   end
end
Prob.probType = checkType('exp');

% MODIFICATION LOG
%
% 981025  hkh  Change SepAlg to struct NLLS