function y = min(a,b,dim)
% tomSym/min - Overloaded

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

% TODO: Fix a nicer version of min, by copying/modifying max.

if nargin==1
    y = -max(-a);
elseif nargin==2
    y = -max(-a,-b);
else %nargin==3
    y = -max(-a,b,dim);
end
