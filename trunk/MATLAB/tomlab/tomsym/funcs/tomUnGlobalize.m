function d = tomUnGlobalize(name,i,id) %#ok
% tomUnGlobalize - Recover data from global variable
%
% d = tomUnGlobalize(d), where d is a tomSym object, reverses the effect
% of tomGlobalize on d, by copying all global data referenced by d into d
% itself.
%
% d = tomUnGlobalize(name,i,id), where i is a positive integer, recovers
% global data set number i (ignoring the "name" input).
%
% See also: tomGlobalize

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-12-23 by rutquist for TOMLAB release 7.7

global tomSymDataStore

if tomSymDataStore{1}~=id
    error('tomSymDataStore no longer contains the requested data.');
end

d = tomSymDataStore{i};
