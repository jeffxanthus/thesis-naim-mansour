function o = rdivide(num,div)
% tomCmplx/rdivide - Division of complex-valued objects

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

a = real(num);
b = imag(num);
c = real(div);
d = imag(div);

q = c.*c + d.*d;

if isnumeric(d) && nnz(d)==0
    o = tomCmplx(a./c, b./c);
else
    o = tomCmplx((a.*c+b.*d)./q, (b.*c - a.*d)./q);
end

