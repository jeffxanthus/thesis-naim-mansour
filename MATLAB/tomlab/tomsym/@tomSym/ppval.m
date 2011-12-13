function y = ppval(pp,xx)
% tomSym/ppval - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if isstruct(xx) && isa(pp,'tomSym')
    warning('tomSymPpval:reversearguments','Argments of ppval call have been reversed.');
    [pp,xx] = deal(xx,pp);
end

if length(pp.dim)==1
    if pp.dim==1
        m = size(xx,1);
        n = size(xx,2);
    else
        if ~any(size(xx)==1)
            error('TomSym currently only handles matrices of dimension <= 2.');
        end
        m = pp.dim;
        n = length(xx);
        if isfield(pp,'orient') && strcmp(pp.orient,'first')
            [m, n] = deal(n,m);
        end
    end
else
    error('TomSym currently only handles matrices of dimension <= 2.');
end

% Simplifications
if all(pp.coefs(:)==0)
    y = zeros(m,n);
else
    y = tomSym(mfilename,m,n,pp,xx);
end
