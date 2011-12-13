function pp = sym2pp(x,pts,flag)
% sym2pp - Convert a tomSym object to a piecewise polynomial
%
% pp = sym2pp(x,pts) returns a fourth order piecewise polynomial which
% takes the same value as x in the points pts, and which has a continuous
% first derivative that is also equal to that of x in those points.
%
% x must be a tomSym expression involving exactly one scalar symbol.
%
% pts is typically a vector of linearly spaced values. The shorter the
% distance between the points is, the higher the accuracy of the
% interpolation will be.
%
% Note that pp is not a tomSym object. It can be evaluated using ppval.
% The points pts will become the breaks in pp.
%
% y = sym2pp(x,pts,'sym') returns a tomSym object. This is useful for
% converting large, slow tomsym expressions into faster interpolation
% tables.
%
% NOTE: sym2pp has low order of accuracy. To obtain an exact representation
% of a tomSym expression for use in other Matlab code, use mcode or mfile.
%
% TIP: sym2pp followed by ppderivative(pp,-1) can be a convenient way to 
% compute an approximate primitive function (antiderivative).
%
% See also: ppval, interp1, mcode, mfile, ppderivative

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

pts = unique(pts);

s = symbols(x,'vector');

if length(s) ~= 1
    error('x must be an expression of a single, scalar symbol.');
end

dx = derivative(x,s);

M = zeros(2,length(pts));
for i=1:length(pts);
    M(1,i) = subs(x,s,pts(i));
    M(2,i) = subs(dx,s,pts(i));
end

no = find(~isfinite(sum(M)));

if ~isempty(no)
    warning('tomsym:infnanpoint',...
        'Expression (or derivative) evaluated to Inf or NaN in some point(s).');
    M(:,no) = [];
    pts(no) = [];
end

if isempty(pts)
    error('All points resulted in Inf or NaN.');
end

dp = diff(pts);

for i=1:length(pts)-1
    d = dp(i);
    B = [1 d; 0 1];
    A = [d^2 d^3; 2*d 3*d^2];
    M(3:4,i) = A\(M(1:2,i+1)-B*M(1:2,i));
end

coefs = M(4:-1:1,1:end-1)';

pp = mkpp(pts,coefs);

if nargin>=3 && flag(1)=='s'
    pp = ppval(pp,s);
end
