function p = tomCmplx(re, im)
% tomCmplx/tomCmplx - Class constructor
%
% c = tomCmplx(re,im) creates a complex-valued symbol from two real-valued
% components.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2010 by Tomlab Optimization Inc.
% Last modified 2010-10-27 by rutquist for TOMLAB release 7.7

% The structure of a tomCmplx object is as follows:
%  p.re - the real component
%  p.im - the imaginary component

if nargin==0
    % Return empty matrix as a tomSym
    % (Sometimes needed for Matlab to get a test object.)
    re = [];
    im = [];
end

if nargin==1
    if isa(re,'tomCmplx')
        p = re;
        return
    else
        im = imag(re);
        re = real(re);
    end
end

% Size check
if ~all(size(re)==size(im));
    error('Re/Im size mismatch when creating tomCmplx object');
end

% Type check
if isa(re,'tomCmplx') || isa(im,'tomCmplx')
    error('Wrong type of arguments to tomCmplx')
end

% Establish priority over tomSym
superiorto('tomSym');

% Create a function
p = struct;
p.re = re;
p.im = im;
p = class(p,'tomCmplx');
