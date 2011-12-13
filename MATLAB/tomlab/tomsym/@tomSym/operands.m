function y = operands(o)
% tomSym/operands - Get all operands from a tomSym
%
% y = operands(o) gets all operands from tomSym object o, as a cell array.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

y = cell(1,nOperands(o));
for i=1:length(y)
    y{i} = operand(i,o);
end
