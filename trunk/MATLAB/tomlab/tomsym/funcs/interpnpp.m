function PPN = interpnpp(method, y, varargin)
% INTERPNPP - Generate PPN struct for ppnval
%
% PPN = INTERPN(METHOD, Y, X1, X2, ...)
%
% METHOD must be one of the methods supported by interp1, i.e. 'linear',
% 'cubic', 'spline', etc. The recommended method is 'spline'. 
%
% The 'cubic' interpolation method is the same as 'pchip', and may result
% in non-smooth surfaces, but in return it preservs maxima and minima.
%
% Y must be a matrix of dimension M by N by ..., where M = length(X1), N =
% length(X2), etc.
%
% This function is used by TomSym's overloaded interp2 and interpn
% functions.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010 by Tomlab Optimization Inc.
% Last modified 2010-01-28 by rutquist for TOMLAB release 7.7

PPN = struct;
N = length(varargin);
PPN.breaks = cell(1,N);

sz = size(y);
for i=1:length(sz);
    if length(varargin{i})~=sz(i)
        error('Dimensions mismatch.');
    end
end

pp = interp1(varargin{1}, y(:,:), method, 'pp');

if ~isstruct(pp)
    error('This Matlab version does not seem to support pp.');
end

c = pp.coefs;
szc = zeros(size(sz));
szc(1) = size(c,2);
PPN.breaks{1} = pp.breaks;
szm = zeros(size(sz));
szm(1) = length(pp.breaks);

for i=2:N
    pp = interp1(varargin{i}, reshape(c,sz(i),numel(c)/sz(i)), method, 'pp');
    c = pp.coefs;
    szc(i) = size(c,2);
    PPN.breaks{i} = pp.breaks;
    szm(i) = length(pp.breaks);
end

PPN.coefs = permute(reshape(c,reshape([szm-1;szc],1,2*N)),[1:2:2*N-1, 2:2:2*N]);
