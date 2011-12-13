function [code, tempD] = mcodestr(f)
% tomSym/mcodestr - Convert tomSym to m-code
%
% [code , tempD] = mcodestr(p) converts tomSym p into a character string 
% of the m-code representing p.
%
% All constants that are used in the code are returned in tempD. 
% The returned code may refer to tempD.
%
% This function is mainly intended for display purpouses, as mcode
% generally produces more efficient code.
%
% See also: mcode

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-08 by rutquist for TOMLAB release 7.7

sp = cell(1,length(f.s));

for i=1:length(sp)
    subp = cell(1,length(f.s(i).a));
    for k=1:length(subp)
        ix = f.s(i).a(k);
        if ix>0;
            subp{k} = sp{ix};
        elseif isnumeric(f.d{-ix}) && numel(f.d{-ix}) < 30
            if issparse(f.d{-ix}) && nnz(f.d{-ix})>=0.3*numel(f.d{-ix})
                subp{k} = mat2str(full(f.d{-ix}));
           else
               subp{k} = mat2str(f.d{-ix});
            end
        else
            subp{k} = ['tempD{' num2str(-ix) '}'];
        end
    end
    sp{i} = getMcode(i,f,subp);
end

code = sp{end};
tempD = f.d;
