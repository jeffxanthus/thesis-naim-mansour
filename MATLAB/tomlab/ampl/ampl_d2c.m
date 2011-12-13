% function d2c = ampl_d2c(x, Prob)
%
% TOMLAB gateway routine for the computation of the constraint
% Jacobian matrix
%
% ampl_dc calls the MEX routine spamfunc (sparse) or amplfunc (dense)
% depending on the value of Prob.AMPL.sparseFlag
%
% The global counter variable n_d2c is incremented

% Marcus Edvall, Tomlab Optimization Inc, E-mail: medvall@tomopt.com
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Jul 31, 2003.   Last modified Jan 17, 2005.

function d2c = ampl_d2c(x, lam, Prob)

global n_d2c

v   = zeros(Prob.AMPL.n_con, 1);
lam = [lam;zeros(Prob.AMPL.n_con-length(lam),1)];

if(Prob.AMPL.sparseFlag == 1)
    d2c = spamfunc(lam);
    d2c = -(d2c - spamfunc(v));
else
    d2c = amplfunc(lam);
    d2c = -(d2c - amplfunc(v));
end

% MODIFICATION LOG
%
% 030731 medv Wrote file
% 030731 ango Added n_d2c handling
% 050117 med  mlint review