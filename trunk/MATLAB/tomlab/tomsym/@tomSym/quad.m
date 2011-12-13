function y = quad(fun,a,b,x)
% tomSym/quad - Numeric quadrature
%
% y = quad(fun,a,b,x) - Computes the integral from a to b of fun with
% respect to x.
%
% If fun contains more than one symbol, then the integral is not computed
% until all the remaining symbols have been substituted by constants.
%
% y = quad(fun,a,b) where fun only contains one symbol, computes the
% quadrature with respect to that symbol.
%
% y = quad(fun,a,b) where fun is a function handle , calls that function 
% with symbolic input to obtain a tomSym function. (This syntax is
% compatible with the Matlab function that is overloaded.)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if ~isa(fun,'tomSym') && nargin<4
    x = tom([],1,1);
    fun = feval(fun,x);
end

s  = symbols(fun,'struct');
ls = sort(fieldnames(s));

if ~exist('x','var')
    if length(ls)>1
        error('Integration varable must be specified.');
    else
        x = s.(ls{1});
    end
end

if numel(x)~=1 || ~tomCmp(x,'tom');
    error('Integration must be with respect to a scalar symbol');
end

if isfield(s,operand(1,x))
    s = rmfield(s,operand(1,x));
    ls = sort(fieldnames(s));
else
    y = fun*(b-a);
    return
end

cs = cell(1,length(ls));
for i=1:length(ls);
    cs{i} = s.(ls{i});
end

fs.fun   = fun;
fs.ivar  = x;
fs.svars = cs;

y = quadsubs(fs,a,b,cs{:});
