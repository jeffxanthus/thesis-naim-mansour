function y = tomreplacement(x)
% tomreplacement - Create a replacement for a tomSym symbol
%
% y = tomreplacement(x) creates a symbol which has the same size and
% integer properties as x, but a different name. This is useful if the same
% equation needs to be duplicated and modified so that the solutions will
% be different.
%
% The "subs" function calls tomreplacement when it encounters subjectTo
% statements. 

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

persistent rcounter

if isempty(rcounter)
    rcounter = 0;
end

if ~tomCmp(x,'tom')
    error('tomreplacement only works directly on tomSym symbols.');
end

rcounter = rcounter + 1;

lbl = operand(1,x);
iflag = operand(2,x);

rix = find(lbl=='_');
if ~isempty(rix)
    rix = rix(end);
    if length(lbl) > rix+1 && lbl(rix+1)=='r' && ...
            all(lbl(rix+2:end) >= '0') && all(lbl(rix+2:end) <= '9')
        lbl = lbl(1:rix-1);
    end
end

lbl = [lbl '_r' num2str(rcounter)];

y = tomSym('tom',size(x,1),size(x,2),lbl,iflag);
