function y=horzcat(varargin)
% tomSym/horzcat - Overload the [] operator

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-10-18 by rutquist for TOMLAB release 7.7

% Remove any empty matrices
sz1 = 0;
sz2 = 0;
nonempty = false(1,length(varargin));
for i=1:length(varargin)
    if isempty(varargin{i})
        if sz1==0
            sz1 = size(varargin{i},1);
        elseif sz1~=size(varargin{i},1)
                warning('tomSym:emptyconcatdim',...
                    'Concatenation involves an empty array with an incorrect number of rows.');
        end
    else
        nonempty(i) = true;
        if sz1==0
            sz1 = size(varargin{i},1);
        elseif size(varargin{i},1)~=sz1
            if any(nonempty(1:i-1))
                error('CAT arguments dimensions are not consistent.');
            else
                sz1 = size(varargin{i},1);
                warning('tomSym:emptyconcatdim',...
                    'Concatenation involves an empty array with an incorrect number of rows.');
            end                
        end
        sz2 = sz2 + size(varargin{i},2);
    end
end

v = varargin(nonempty);
i = 1;
while i<=length(v)
    if tomCmp(v{i},'horzcat');
        v = [v(1:i-1), operands(v{i}), v(i+1:end)];
    end
    i = i+1;
end

if length(v)==1
    y = v{1};
else
    y = tomSym(mfilename, sz1, sz2, v{:});
end
