% function [f,g] = smoothmin(x, y, delta)
%
% smoothmin is a smooth approximation to min(x,y) when delta > 0.
% Small delta gives a better approximation (e.g. 1e-3 or 1e-4).
% If abs(x-y) >= delta, there is no error.
% If x = y, the error in f is delta/2.
%
% INPUT:
% x           Value 1
% y           Value 2
% delta       Tolerance
%
% OUTPUT:
% f           Function value
% g           Gradient vector

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written May 24, 2004.  Last modified May 24, 2004.

function [f,g] = smoothmin(x, y, delta)

t = (x-y)/2;

if t > delta/2
    f =  y;
    g = [0; 1];
elseif t < -delta/2
    f =  x;
    g = [1; 0];
else
    f = t^2/delta + delta/4 + (x+y)/2;
    g = [ t/delta + 0.5
        -t/delta + 0.5];
end

% MODIFICATION LOG
%
% 040524  med  Written, Gill's code