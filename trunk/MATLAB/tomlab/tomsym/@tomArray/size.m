function varargout = size(o,i)
% tomArray/size - Get the size of a tomArray object.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

varargout = cell(1,nargout);

if nargout>1
    for i = 1:nargout;
        if i<=length(o.sz)
            varargout{i} = o.sz(i);
        else
            varargout{i} = 1;
        end
    end
    if nargout<length(o.sz)
        varargout{end} = prod(o.sz(nargout:end));
    end
else
    if nargin>1
        varargout{1} = o.sz(i);
    else
        % Stupid Matlab gets confused if length(size(x)) < 2.
        varargout{1} = [o.sz ones(1,max(0,2-length(o.sz)))];        
    end
end
