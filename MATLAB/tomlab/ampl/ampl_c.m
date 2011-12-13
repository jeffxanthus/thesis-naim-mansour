% function c = ampl_c(x, Prob)
%
% TOMLAB gateway routine for the computation of nonlinear constraint values c(x)
%
% ampl_g calls the MEX routine spamfunc (sparse) or amplfunc (dense)
% depending on the value of Prob.AMPL.sparseFlag
%
% The global counter variable n_c is incremented

% Marcus Edvall, Tomlab Optimization Inc, E-mail: medvall@tomopt.com
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Jul 31, 2003.   Last modified Jan 17, 2005.

function c = ampl_c(x, Prob)

global n_c

if(Prob.AMPL.sparseFlag == 1)
    c = spamfunc(x, 3);
else
    c = amplfunc(x, 3);
end

n_c = n_c + 1;

% MODIFICATION LOG
%
% 030731 medv Wrote file
% 030731 ango Added n_c handling
% 050117 med  mlint review