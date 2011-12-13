function o = tomArrayIdx(s,i)
% tomArrayIdx class constructor
%
% tomArrayIdx(s,i) creates a tomArrayIdx with the symbolic name s
% containing the indexes i.
%
% A tomArrary idx is used for indexing the tomArray class, and allows a
% more powerful syntax as compared to when using string constants.
%
% A tomArrayIdx is often used as the mathematical notation "for all"
% For example, the code:
%   k = tomArrayIdx('k',1:3); 
%   c = ( x(k) == y(k+1) )
% can be interpreted as "the constraint c says that x(k) equals y(k+1) for
% all k in (1:3).
%
% See also: tomArray

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-26 by rutquist for TOMLAB release 7.7

persistent autoix

if isempty(autoix)
    autoix = 1;
end

if nargin==0
    o.s = ['ix' num2str(autoix)];
    o.i = ':';
    autoix = autoix+1;
    o = class(o,'tomArrayIdx');
    return
end

if nargin==1
    if isa(s,'tomArrayIdx')
        o = s;
    elseif isvarname(s)
        o.s = s;
        o.i = ':';
        o = class(o,'tomArrayIdx');
    elseif isnumeric(s)
        o.s = ['ix' num2str(autoix)];
        o.i = s;
        autoix = autoix+1;
        o = class(o,'tomArrayIdx');
    else
        error([class(s) 'cannot be converted into a tomArrayIdx.']);
    end
else
    if ~(isvarname(s) && isnumeric(i))
        error('Wrong type input argument to tomArrayIdx.');
    end
    o.s = s;
    o.i = i;
    o = class(o,'tomArrayIdx');
end
