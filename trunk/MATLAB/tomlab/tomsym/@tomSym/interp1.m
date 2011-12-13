function yi = interp1(x,y,xi,method,extrap) %#ok
% tomSym/interp1 - Overloaded function

% The "extrap" argument is not used, but acts as a placeholder in case the
% user would supply it. (Extrapolation is the default in any case.)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-08-13 by rutquist for TOMLAB release 7.7

if nargin<4
    method = 'linear';
end

if nargin==2 || (nargin==3 && ischar(xi))
    if nargin==3
        method = xi;
    end
    xi = y;
    y  = x;
    if size(y,2)>1
        x = 1:size(y,2);
    else
        x = 1:length(y);
    end
end

if length(x) <= 1
    error('tomSym:interp1:NotEnoughPts', 'There should be at least two data points.')
end

if ~((length(x) == length(y) && length(y)==numel(y)) || length(x) == size(y,1))
    error('tomSym:interp1:Mismatch', 'x and y must have the same length.')
end
    
if ( tomCmp(x,'smtimes') && (tomCmp(xi,'smtimes') || ...
        (tomCmp(xi,'mtimes') && numel(xi)==1)) || ...
        tomCmp(x,'plus') && tomCmp(xi,'plus') ) && ...
        isequal(operand(1,x),operand(1,xi))
    % Simplifications necessary for atPoints to work as expected.
    yi = interp1(operand(2,x),y,operand(2,xi),method);
    return;
end

if isnumeric(xi) && numel(xi)==1 && isnumeric(x)
    % Simplify for exact match.
    idx = find(x==xi);
    if ~isempty(idx)
        yi = lookup(y,idx);
        return
    end
end

if isequal(vec(x),vec(xi))
    yi = reshape(y,size(xi));
elseif isa(x,'tomSym') || isa(y,'tomSym')
    if strcmpi(method,'linear')
        yi = interp1l(x,y,xi);
    elseif strcmpi(method,'spline');
        yi = interp1s(x,y,xi);
    else
        error('This kind of interpolation is currently not supported by tomSym.');
    end
else
    switch lower(method)
        case 'linear'
            yi = interp1l(x,y,xi);
        case 'spline'
            yi = interp1s(x,y,xi);
        otherwise
            pp = interp1(x,y,method,'pp');
            if ~isstruct(pp);
                % TODO: Fix this to be compatible with Matlab 6.5
                error('tomSym:interp1:OldMatlabVersion',...
                    'This method of interpolation is not supported by tomSym under Matlab 6.5 or earlier');
            end
            yi = ppval(pp,xi);
    end
end
