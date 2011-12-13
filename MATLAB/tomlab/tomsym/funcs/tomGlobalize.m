function d = tomGlobalize(d,name)
% tomGlobalize - Store data in global variable, for memory efficency
%
% NOTE: tomGlobalize is a temporary fix for certain memory-issues in 
% older versions of Matlab. It may be removed in a future TOMLAB release.
% The recommended way to make a symbolic reference out of a constant 
% is: d = tomSym(d)
%
% d = tomGlobalize(d) turns d into a tomData symbol, and stores the actual 
% data associated with d in the global variable tomSymDataStore.
%
% The created symbol uses very little memory, which will reduce memory
% requirements for operations that involve multiple copies of the data.
% If a symbol contains, for example, both d and sin(d), then only d
% will be stored in memory, and sin(d) evaluated at run-time. (If d is not
% globalized, then sin(d) will be treated as a constant and stored in
% memory with the symbol.)
%
% d = tomGlobalize(d,'name') assigns a display name to the symbol. (The 
% default is the name of the input variable.)
%
% tomGlobalize('clear') removes all previously stored global data, freeing
% up memory.
%
% The effect of tomGlobalize can be undone by applying tomUnGlobalize or
% eval.
%
% Note: If Symoblic objects that involve tomData are saved using
% Matlab's "save" command, the data will not be included, and subsequent
% loads will yield inconsistent objects. Always use tomUnGlobalize on 
% symbolic objects before they are saved.
%
% See also: tomUnGlobalize

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-12-23 by rutquist for TOMLAB release 7.7

global tomSymDataStore

if isempty(tomSymDataStore)
    tomSymDataStore = {rand(1)};
end

if ischar(d) && strcmpi(d,'clear')
    tomSymDataStore = {};
    return
end

if nargin<2
    name = inputname(1);
end

for i=1:length(tomSymDataStore)
    if isequalwithequalnans(d,tomSymDataStore{i})
        if isempty(name)
            name = ['tempG' num2str(i)];
        end
        d = tomSym('tomUnGlobalize',size(d,1),size(d,2),name,i,tomSymDataStore{1});
        return
    end
end

tomSymDataStore{end+1} = d;
i = length(tomSymDataStore);

if isempty(name)
    name = ['tempG' num2str(i)];
end

d = tomSym('tomUnGlobalize',size(d,1),size(d,2),name,i,tomSymDataStore{1});
