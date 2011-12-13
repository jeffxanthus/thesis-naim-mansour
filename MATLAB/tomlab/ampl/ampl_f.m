% function f = ampl_f(x, Prob)
%
% TOMLAB gateway routine for the computation of function values f(x)
%
% ampl_f calls the MEX routine spamfunc (sparse) or amplfunc (dense)
% depending on the value of Prob.AMPL.sparseFlag
%
% The global counter variable n_f is incremented

% Marcus Edvall, Tomlab Optimization Inc, E-mail: medvall@tomopt.com
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Jul 31, 2003.   Last modified Jan 17, 2005.

function f = ampl_f(x, Prob)

global n_f

if(Prob.AMPL.sparseFlag == 1)
    f = spamfunc(x, 2);
else
    f = amplfunc(x, 2);
end

f = Prob.AMPL.objType * f;

n_f = n_f + 1;

% MODIFICATION LOG:
%
% 030731 med  Wrote file
% 030731 ango Added n_f handling
% 050117 med  mlint review