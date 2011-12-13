% function f = ego_f(x, Prob)
%
% Expected improvement function

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2006 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Sep 29, 1998.   Last modified August 24, 2009.

function f = ego_f(x, Prob)

f = cgolib(204, Prob.EGO.CGOLIB.daceid(Prob.EGO.CGOLIB.TRANSFORM+1), x);

switch(Prob.CGO.EITRANSFORM)
    case 0,
       % f = f; % f is already defined
    case 1,
       if f>0 
          f = 1e10;
       else
          f = -log(-f+1e-5);
       end 
    case 2,
        if f>0
           f = 1e10;
        else
            % Is 1e-5 an apropriate margin?
            f = -1/(f+1e-5);
        end
end

% MODIFICATIONS
% 040930  hkh  Safe guard, set f=1E10 if not finite
% 050218  hkh  Clean up
% 051112  frhe Everything removed, and cgolib_mex call added.
% 060105  frhe Added -log(-f) transformation, and safe for values > 0. 
% 060207  ango Change cgolib_mex -> cgolib
% 090824  hkh  Avoid f=f;
