function [m,n] = size(o,i)
% tomSym/size - Get the size of a tomSym object

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if nargout>1
    m = o.s(end).sz1;
    n = o.s(end).sz2;
else
    if nargin>1
        if(i==1)
            m = o.s(end).sz1;
        elseif i==2
            m = o.s(end).sz2;
        else
            error('TomSym objects only have two dimensions.');
        end
    else
        m = [o.s(end).sz1 o.s(end).sz2];
    end
end
