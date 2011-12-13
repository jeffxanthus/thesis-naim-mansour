function A = subsasgn(A, S, B)
% tomSym/subsasgn - Subscripted assignemnt, overloaded for tomSym.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-03-04 by rutquist for TOMLAB release 7.7

if ~strcmp(S.type,'()')
    error('Invalid type of assignment for tomSym');
end

% Convert to double if possible
if tomCmp(A,'constant')
    A = double(A);
end
if tomCmp(B,'constant')
    B = double(B);
end

if ~isnumeric(S.subs{1})
    if ischar(S.subs{1}) && strcmp(S.subs{1},':')
        % Expand : to vector of indexes
        if isempty(A)
            S.subs{1} = 1:size(B,1);
        else
            S.subs{1} = 1:size(A,1);
        end
    elseif islogical(S.subs{1})
        S.subs{1} = find(S.subs{1});
    else
        error(['Invalid index for tomSym: ' class(S.subs{1})]);
    end
end

if length(S.subs)>1 && ~isnumeric(S.subs{2})
    if ischar(S.subs{2}) && strcmp(S.subs{2},':')
        if isempty(A)
            S.subs{2} = 1:size(B,2);
        else
            S.subs{2} = 1:size(A,2);
        end
    else
        error(['Invalid index for tomSym: ' class(S.subs{2})]);
    end
end

if length(S.subs)>2
        error('Too many indexes for tomSym.');
end

% Determine size of result
if length(S.subs)==1
    if size(A,2)==1 && size(A,1)>1
        sz2 = 1;
        sz1 = max(S.subs{1});
    elseif size(A,1)==1
        sz1 = 1;
        sz2 = max(S.subs{1});
    else
        if max(S.subs{1}) > numel(A)
            error('In an assignment  A(I) = B, the matrix A cannot be resized.');
        end
        sz1 = size(A,1);
        sz2 = size(A,2);
    end
else %length(S.subs)==2
    sz1 = max(S.subs{1});
    sz2 = max(S.subs{2});
end
sz1 = max(sz1,size(A,1));
sz2 = max(sz2,size(A,2));

if length(S.subs)==1
    if size(A,2)==1 && size(A,1)>1
        jj = ones(size(S.subs{1}));
        ii = S.subs{1};
    elseif size(A,1)==1
        ii = ones(size(S.subs{1}));
        jj = S.subs{1};
    else
        [ii,jj] = ind2sub([sz1,sz2],S.subs{1});
    end
else
    [jj,ii] = meshgrid(S.subs{2},S.subs{1});
end

if numel(ii) ~= numel(B)
    if numel(B)==1
        B = repmat(B,size(ii));
    else
        error('Subscripted assignment dimension mismatch.');
    end
end

if isempty(A)
    if length(S.subs)==1 && all(S.subs{1}==1:numel(B))
        A = vec(B)';
    elseif sz1==length(S.subs{1}) && all(S.subs{1}==1:sz1) && ...
            sz2==length(S.subs{2}) && all(S.subs{2}==1:sz2)
        A = B;
    else
        A = sparse(ii,jj,B,sz1,sz2);
    end
else
    if sz1>size(A,1) || sz2>size(A,2)
        if tomCmp(A,'sparse')
            A = sparse(operand(1,A),operand(2,A),operand(3,A),sz1,sz2);
        elseif tomCmp(A,'setdiag')
            A = sparse(1:size(A,1),1:size(A,1),operand(1,A),sz1,sz2);
        else
            [ja,ia] = meshgrid(1:size(A,2),1:size(A,1));
            A = sparse(ia,ja,A,sz1,sz2);
        end
    end
    if ~any(lookup(pattern(A),sub2ind([sz1,sz2],ii,jj)))
        A = A + sparse(ii,jj,B,sz1,sz2);
    else
        if tomCmp(A,'tom')
            disp('Tip: To create a symbolic object that is a mixture of unknowns and constants,');
            disp('use the concatenation operators. For example: x = [0; tom(''x2_to_x4'', 3, 1)]');            
        end
        error('Changing non-zero elements of a tomSym array is not possible.');
    end
end
