% function z = nlp_d2Lv(x,lam,v,Prob)
%
% Callback function used by KNITRO. Calculates the Hessian of the 
% Lagrangian (f(x)+lam*c(x)) and multiplies the result with the vector v.

% Anders Goran, Tomlab Optimization Inc, E-mail: anders@tomopt.com
% Copyright (c) 2003 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written May 13, 2003.  Last modified Jun 10, 2003.

function z = nlp_d2Lv(x,lam,v,Prob)
z = nlp_d2L(x,lam,2,Prob) * v;

% MODIFICATION LOG
%
% 030513 ango Wrote file
% 030610 ango Revise comments