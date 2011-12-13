function mad_derivs=getinternalderivs(mad_obj)
%
%getinternalderivs   
%       Get sparse derivatives from fmad object.
%       fmad is Automatic Differentiation with forward mode.
%
%       See MAD users manual for detail.
%
%       ALSO SEE  fmad, getderivs
%
% 

nargin;
disp('DUMMY for TOMLAB if MAD not installed');
mad_derivs=[];