function y=setrows(n,varargin)
% setrows - Overloaded function
%
% y = ekron(n,ix1,M1,ix2,M2,...) sets the rows given by ix1 to M1, and the
% rows given by ix2 to M2, etc. This is equivalent to:
%
% y = sparse(n,size(M1,2));   
% y(ix1,:) = y(ix1,:) + M1;
% y(ix2,:) = y(ix2,:) + M2;
% ...
%
% This is a specialized form of sparse matrices with only a few nonzero
% rows. 
%
% TomSym uses setrows for reasons of memory efficiency.
%
% See also: submatrix

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2011 by Tomlab Optimization Inc.
% $Id$

if length(varargin) < 2
    error('Too few arguments to setrows');
end

if mod(length(varargin),2) ~= 0
    error('Wrong number of arguments to setrows');
end

if length(varargin)==2 && length(varargin{1})==n
    % One full matrix.
    inverseix(varargin{1}) = 1:length(varargin{1});
    y = submatrix(varargin{2},inverseix,':');
    return
end

sz2 = size(varargin{2},2);

usedby = zeros(1,n);

for i=1:length(varargin)/2
    ix = varargin{2*i-1};
    if length(ix)~=size(varargin{2*i},1)
        error('Matrix size must match index list length.');
    end
    xi = find(usedby(ix));
    ui = usedby(ix(xi));
    usedby(ix) = i;
    for k=1:i-1
        if any(ui==k);
            xi1 = usedby(varargin{2*k-1})==i;
            xi2 = find(ui==k);
            varargin{2*i} = varargin{2*i} + ...
                setrows(size(varargin{2*i},1),xi2,submatrix(varargin{2*k},find(xi1),':'));
            varargin{2*k-1}(xi1) = [];
            varargin{2*k} = submatrix(varargin{2*k},find(~xi1),':');
        end    
    end
end

isok = true(1,length(varargin)/2);
for i=1:length(varargin)/2
    if isempty(varargin{2*i}) || iszero(varargin{2*i})
        isok(i) = false;
    end
end

isok = reshape([isok; isok],1,length(varargin));
y = tomSym(mfilename,n,sz2,n,varargin{isok});
