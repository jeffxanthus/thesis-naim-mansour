function y = fliplr(a)
% tomSym/fliplr - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

% Simplifications
if size(a,2) == 1
    y = a;
elseif strcmp(operator(a),'fliplr')
    y = operand(1,a);
elseif strcmp(operator(a),'horzcat')
    n = nOperands(a);
    fa = cell(1,n);
    for i=n:-1:1
      fa{n-i+1} = fliplr(operand(i,a));
    end
    y = horzcat(fa{:});
else
    y = tomSym(mfilename,size(a,1),size(a,2),a);
end
