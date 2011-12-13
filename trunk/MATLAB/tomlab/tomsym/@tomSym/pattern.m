function p = pattern(f)
% tomSym/pattern - The sparsity pattern of a tomSym object
%
% p = pattern(f) computes the (possibly sparse) logical matrix p, that is
% true (1) wherever f might be nonzero, and false (0) elsewhere. 
%
% The pattern might contain more nonzeros than necessary, but preserves
% sparsity patterns for function that map zero to zero. 

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

sp = cell(1,length(f.s));
dp = cell(1,length(f.d));

for i=1:length(dp)
    if isnumeric(f.d{i})
        dp{i} = f.d{i}~=0;
    elseif islogical(f.d{i})
        dp{i} = f.d{i};
    else
        dp{i} = [];
    end
end

v = lastdepvec(f);

for i=1:length(sp)
    subp = cell(1,length(f.s(i).a));
    for k=1:length(subp)
        ix = f.s(i).a(k);
        if ix>0;
            subp{k} = sp{ix};
        else
            subp{k} = dp{-ix};
        end
    end
    sp{i} = getPattern(i,f,subp);
    clr = find(v==i);
    for j=1:length(clr)
        sp{clr(j)} = []; % Free up some memory
    end
end

p = sp{end};
