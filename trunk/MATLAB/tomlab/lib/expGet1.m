% expGet1 picks up four parameters for exponential fitting from Prob structure
%
% function [SepAlg, p, wType, eType, infCR] = expGet1(Prob);
%
% INPUT:
%  Prob     Problem structure
%
% OUTPUT:
%  SepAlg   Flag if Separable NLLS solution
%  p        No of exponential terms
%  wType    Weighting type
%  eType    Type of exponential terms
%  infCR    Information criteria for selection of best number of terms
               
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Feb 17, 1998.  Last modified Nov 5, 2000.

function [SepAlg, p, wType, eType, infCR] = expGet1(Prob)

if isempty(Prob)
   SepAlg=0;
   p=     2;
   wType= 0;
   eType= 1;
   infCR= 0;
else
   SepAlg=Prob.LS.SepAlg;
   p=     Prob.ExpFit.p;
   wType= Prob.ExpFit.wType;
   eType= Prob.ExpFit.eType;
   infCR= Prob.ExpFit.infCR;
end

% MODIFICATION LOG
%
% 981025  hkh  Change SepAlg to struct NLLS
% 001105  hkh  Change structure NLLS to LS