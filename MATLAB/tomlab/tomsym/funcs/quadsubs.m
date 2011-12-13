function y = quadsubs(fs,a,b,varargin)
% tomSym/quad - Numeric quadrature with symbol substitution
%
% y = quad(fs,a,b,subs1,...) - Computes the integral from a to b of fun,
% substituting values subs1, subs2, ... for the symbols that it contains.
%
% fs must be a struct with the fields:
%     fun   - a tomSym object representing the funtion to be integrated
%     ivar  - a tomSym symbol identifying the integration variable
%     svars - a cell list of tomSym symbols defining the order of substvals
%
% This function is typically not called directly, but calls will be
% generated as a consequence of calls to tomSym/quad.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-10-02 by rutquist for TOMLAB release 7.7

% TODO: Substitute subtrees instead of symbols, when possible. Identify
% minmal cut?

if numel(a)~=1 || numel(b)~=1
    error('Integration boundaries must be scalar-');
end

if isequal(a,b)
    y = zeros(size(fs.fun));
    return
end

if ~isa(fs.fun,'tomSym')     
    y = (b-a)*fs.fun;
    return
end

ss = struct;
for i=1:length(varargin);
    ss.(operand(1,fs.svars{i})) = varargin{i};
end
fs.fun = subs(fs.fun,ss);

if isnumeric(fs.fun) || ~isdependent(fs.fun, fs.ivar)
    y = (b-a)*fs.fun;
    return
end

s  = symbols(fs.fun,'struct');

if isfield(s,operand(1,fs.ivar))
    s = rmfield(s,operand(1,fs.ivar));
    ls = sort(fieldnames(s));
else
    y = fs.fun*(b-a);
    return
end

cs = cell(1,length(ls));
for i=1:length(ls);
    cs{i} = s.(ls{i});
end
fs.svars = cs;

if ~(numel(fs.ivar)==1 && tomCmp(fs.ivar,'tom'));
    error('Integration must be with respect to a scalar symbol');
end

if ~isempty(cs) || isa(a,'tomSym') || isa(b,'tomSym')
    y = tomSym(mfilename,size(fs.fun,1),size(fs.fun,2),fs,a,b,cs{:});
else
    y = quadv(@(xx) subs(fs.fun,fs.ivar,xx),a,b);
end

