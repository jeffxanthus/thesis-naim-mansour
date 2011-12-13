function o = mrdivide(num,div)
% tomCmplx/mrdivide - Matrix division of complex-valued objects

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

if isnumeric(div) && size(div,1) == size(div,2)
    o = num*inv(div);
    return
end

if numel(num)==1 && numel(div)==1
    o = num./div;
    return
end

error('Matrix division with symbolic complex divisor is not supported.');
