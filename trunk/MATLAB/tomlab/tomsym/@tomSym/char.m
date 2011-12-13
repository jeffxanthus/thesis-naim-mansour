function c = char(f)
% tomSym/char - Convert tomSym to char.
%
% c = char(p) converts tomSym p into a character string.
%
% See also: mcodestr

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-19 by rutquist for TOMLAB release 7.7

sp = cell(1,length(f.s));
dp = cell(1,length(f.d));

for i=1:length(dp)
    dp{i} = tochar(f.d{i});
end

for i=1:length(sp)
    if i>length(sp)-30
        subp = cell(1,length(f.s(i).a));
        for k=1:length(subp)
            ix = f.s(i).a(k);
            if ix>0;
                subp{k} = sp{ix};
            else
                subp{k} = dp{-ix};
            end
        end
        sp{i} = getChar(subsymb(i,f),subp);
    else
        sp{i} = ['[... ' num2str(f.s(i).sz1) ' by ' num2str(f.s(i).sz2) ' tomSym ...]'];
    end
end

c = sp{end};
