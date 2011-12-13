% function g = ampl_g(x, Prob)
%
% TOMLAB gateway routine for the computation of gradient values g(x)
%
% ampl_g calls the MEX routine spamfunc (sparse) or amplfunc (dense)
% depending on the value of Prob.AMPL.sparseFlag
%
% The global counter variable n_g is incremented

% Marcus Edvall, Tomlab Optimization Inc, E-mail: medvall@tomopt.com
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Jul 31, 2003.   Last modified Jan 17, 2005.

function g = ampl_g(x, Prob)

global n_g

if(Prob.AMPL.sparseFlag == 1)
    g = spamfunc(x, 4);
else
    g = amplfunc(x, 4);
end

g = Prob.AMPL.objType * g;

n_g = n_g + 1;

% MODIFICATION LOG
%
% 030731 medv Wrote file
% 030731 ango Added n_g handling
% 050117 med  mlint review