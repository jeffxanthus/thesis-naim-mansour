function [c, n] = flattencellarray(c,n0)
% flattencellarray - Flatten cell array-of-arrays.
%
% [c, n] = flattencellarray(c,n0) returns a flattened cell array c and a
% "name-list" n, starting with the base name n0.
%
% The name list will contain references to the originall array-of-arrays
% for each cell. 
%
% For example: if c = {1, {2, 3}} and n0 = 'c' then n{3} will be 'c{2}{2}'.
%
% If c contains tomSym expressions of the type a <= b <= c, then those will
% also be split up and named as 'c{i}<1>' and 'c{i}<2>', etc.


% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

if nargin < 2
    n0 = '';
end

c = c(:)';
n = cell(size(c));

i = 1;
k = 1;
while i <= length(c)
    if iscell(c{i})
        [c1, n1] = flattencellarray(c{i}, [n0, '{', num2str(k), '}']);
        c = [c(1:i-1) c1 c(i+1:end)];
        n = [n(1:i-1) n1 n(i+1:end)];
        i = i+length(n1);
    elseif tomCmp(c{i},'le') && tomCmp(operand(1,c{i}),'le')
        cc = c{i};
        ce = {};
        while tomCmp(cc,'le')
            ce{end+1} = operand(2,cc);
            cc = operand(1,cc);
        end
        ce{end+1} = cc;
        ce = ce(end:-1:1);
        c1 = cell(1,length(ce)-1);
        n1 = cell(1,length(ce)-1);
        for r=1:length(ce)-1;
            c1{r} = (ce{r} <= ce{r+1});
            n1{r} = [n0 '{', num2str(k), '}', '<', num2str(r), '>'];
        end
        c = [c(1:i-1) c1 c(i+1:end)];
        n = [n(1:i-1) n1 n(i+1:end)];
        i = i+length(n1);
    else
        n{i} = [n0, '{', num2str(k), '}'];
        i = i+1;
    end
    k = k+1;
end
