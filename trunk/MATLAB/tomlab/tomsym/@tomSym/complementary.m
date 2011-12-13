function y = complementary(a,b)
% complementary - overloaded function
%
% complementary(a,b) creates a complementary codition:
%
% For each element i in a and b, either a(i)==0 or b(i)==0.
%
% Complementary conditions require the KNITRO solver.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if ~all(size(a)==size(b))
    error('Complementary condition requires that arrays have the same size.');
end

%if isnumeric(a)
%    fi = find(a);
%    if isempty(fi)
%        y = [];
%    else
%        y = ( b(fi)==0 );
%    end
%elseif isnumeric(b)
%    fi = find(b);
%    if isempty(fi)
%        y = [];
%    else
%        y = ( a(fi)==0 );
%    end
%else
    y = tomSym(mfilename,size(a,1),size(a,2),a,b);
%end
