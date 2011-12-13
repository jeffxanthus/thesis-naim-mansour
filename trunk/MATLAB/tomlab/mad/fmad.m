function mad_x=fmad(x,dx)
%
%fmad   Create fmad object for calculcation of function and derivatives
%       It computes sparse matrcies of a general scalar nonlinear mapping.        
%       fmad is Automatic Differentiation with forward mode.
%
%       See MAD users manual for detail.
%
%       ALSO SEE  getinternalderivs, getderivs
%
% 

nargin;
disp('DUMMY for TOMLAB if MAD not installed');
mad_x=[];