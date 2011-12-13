function y=setrows(n,varargin)
% setrows - Set rows in an otherwise zero matrix.
%
% y = ekron(n,ix1,M1,ix2,M2,...) sets the rows given by ix1 to M1, and the
% rows given by ix2 to M2, etc. This is equivalent to:
%
% y = sparse(n,size(M1,2));   
% y(ix1,:) = y(ix1,:) + M1;
% y(ix2,:) = y(ix2,:) + M2;
% ...
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

used = false(1,n);
ok = true;
for i=1:length(varargin)/2
    if any(used(varargin{2*i-1}))
        ok = false;
        break;
    end
end

if ~issparse(varargin{2})
    y = zeros(n,size(varargin{2},2));
    if ok
        for i=1:length(varargin)/2
            y(varargin{2*i-1},:) = varargin{2*i};
        end
    else
        for i=1:length(varargin)/2
            y(varargin{2*i-1},:) = y(varargin{2*i-1},:) + varargin{2*i};
        end
    end
else
    if ok
        %ic = cell(1,length(varargin)/2);
        %jc = cell(1,length(varargin)/2);
        %vc = cell(1,length(varargin)/2);
        %for i=1:length(varargin)/2
        %    [i1,j1,v1] = find(varargin{2*i});
        %    ic{i} = varargin{2*i-1}(i1);
        %    jc{i} = j1;
        %    vc{i} = v1;
        %end
        %y = sparse(vertcat(ic{:}),vertcat(jc{:}),vertcat(vc{:}),n,size(varargin{2},2)); 
        y = sparse(n,size(varargin{2},2));
        for i=1:length(varargin)/2
            y(varargin{2*i-1},:) = varargin{2*i};
        end
    else
        y = sparse(n,size(varargin{2},2));
        for i=1:length(varargin)/2
            y(varargin{2*i-1},:) = y(varargin{2*i-1},:) + varargin{2*i};
        end
    end
end
