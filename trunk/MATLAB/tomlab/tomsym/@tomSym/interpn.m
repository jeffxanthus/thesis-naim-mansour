function vi = interpn(varargin)
% TOMSYM/INTERPN - overloaded function
%
% VI = INTERPN(X1,X2,X3,...,V,Y1,Y2,Y3,...) interpolates V to find VI.
%
% VI = INTERP2(X1,X2,X3,...,V,Y1,Y2,Y3,...,METHOD) specifies the method to
% use. 
%
% X1, X2, ... and V must be numeric. Only Y1, Y2, ... are currently
% implemented as tomSym variables.
%
% Extrapolation is currently handled differently as compared to the
% standard Matlab function. The same method is used for both interpolation
% and extrapolation. This may change in the future. It is recommended to
% use a larger table and interpolate, rather than to extrapolate.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010 by Tomlab Optimization Inc.
% Last modified 2010-01-28 by rutquist for TOMLAB release 7.7

if ischar(varargin{end})
    method = varargin{end};
    Ni = length(varargin)-1;
else
    method = 'linear';
    Ni = length(varargin);
end

N = (Ni-1)/2;
if N~=round(N)
    error('Incorrect number of arguments');
end
    
if method(1) == 'c'
    method = 'v5cubic';
end

for i=1:N
    if ~isnumeric(varargin{i})
        error('Symoblic X arguments have not yet been implemented.');
    end
    if ~all(size(varargin{N+1+i})==size(varargin{N+1+1}))
        error('Y-sizes mismatch');
    end
end

if isnumeric(varargin{N+1})
    PPN = interpnpp(method, varargin{N+1}, varargin{1:N});
    vi = tomSym('ppnval', size(varargin{N+1+1},1), ...
        size(varargin{N+1+1},2), PPN, varargin{N+1+1:Ni});
else
    error('Not all variations of interpn have been implemented for tomSym yet.');
end

