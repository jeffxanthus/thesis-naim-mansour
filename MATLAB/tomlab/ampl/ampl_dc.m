% function cjac = ampl_dc(x, Prob)
%
% TOMLAB gateway routine for the computation of the constraint
% Jacobian matrix
%
% ampl_dc calls the MEX routine spamfunc (sparse) or amplfunc (dense)
% depending on the value of Prob.AMPL.sparseFlag
%
% The global counter variable n_dc is incremented

% Marcus Edvall, Tomlab Optimization Inc, E-mail: medvall@tomopt.com
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Jul 31, 2003.   Last modified Jan 17, 2005.

function cjac = ampl_dc(x, Prob)

global n_dc

if(Prob.AMPL.sparseFlag == 1)
    cjac = spamfunc(x, 5);
else
    cjac = amplfunc(x, 5);
end

n_dc = n_dc + 1;

% MODIFICATION LOG
%
% 030731 medv Wrote file
% 030731 ango Added n_dc handling
% 050117 med  mlint review