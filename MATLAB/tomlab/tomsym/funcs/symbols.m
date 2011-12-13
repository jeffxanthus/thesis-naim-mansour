function list = symbols(o,flag)
% symbols - List all symbols used in a tomSym
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
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

% If o is a tomSym then the overloaded function is called. This is just a
% placeholder function which returns an empty list of symbols.

if nargin<2
    flag = 'list';
end

if isa(o,'cell')
    o = tomSym('cellarray',size(o,1),size(o,2),o{:});
    list = symbols(o,flag);
    return
end

switch flag
    case 'struct'
        list = struct();
    case 'vector'
        [];
    case 'all'
        list = struct();
    case 'list'
        list = {};
    otherwise
        error('Unknown flag');
end
