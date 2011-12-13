% function H = ampl_H(x, Prob)
%
% TOMLAB gateway routine for the computation of the Hessian H(x)
%
% ampl_g calls the MEX routine spamfunc (sparse) or amplfunc (dense) 
% depending on the value of Prob.AMPL.sparseFlag

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Jul 31, 2003.   Last modified Jan 17, 2005.

function W = ampl_H(x, Prob)

% The Lagrangian multipliers are inputs to the Hessian calculation.
v = zeros(Prob.AMPL.n_con,1);
if(Prob.AMPL.sparseFlag == 1)
    W = spamfunc(v);
else
    W = amplfunc(v);
end

W = Prob.AMPL.objType * W;

% MODIFICATION LOG
%
% 030731 med  Wrote file
% 030731 ango Added n_H handling
% 040701 med  objtype added for W
% 050117 med  mlint review