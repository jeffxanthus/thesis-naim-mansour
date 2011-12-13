function list = symbols(o,flag)
% tomSym/symbols - List all symbols used in a tomSym
%
% list = symbols(o) returns a cell array of symbol names used in the tomSym
% object or cell array o.
%
% s = symbols(o,'struct') returns a struct
% v = symbols(o,'vector') returns a tomSym vector
%
% Only symbols that are needed for evaluation are returned. Symbols whose
% value are given implicitly by subjectTo calls, and indexes used in 
% subsVec are not returned. s = symbols(o,'all') returns a struct including
% these symbols.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if nargin<2
    flag = 'list';
end

switch flag
    case 'struct'
        list = symbstruct(struct, o, false);
    case 'vector'
        l = symbstruct(struct, o, false);
        fl = sort(fieldnames(l));
        c = cell(1, length(fl));
        for i=1:length(fl)
            c{i} = vec(l.(fl{i}));
        end
        list = vertcat(c{:});
    case 'all'
        list = symbstruct(struct, o, true);
    case 'list'
        list = sort(fieldnames(symbstruct(struct, o, false)));
    otherwise
        error('Unknown flag');
end
