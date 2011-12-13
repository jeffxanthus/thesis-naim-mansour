function M = tomArray(X,sz,ni,intflag)
% tomArray class constructor
%
%  M = tomArray(X)       converts X into a tomArray
%  M = tomArray(X,sz)    converts X into tomArray form with size = sz.
%  M = tomArray(X,sz,ni) also asigns index names ni.
%  M = tomArray('x',sz)  creates a tomsym named 'x' as basis for the array.
%  M = tomArray('x',sz,ni,'int') creates an integer tomSym as basis.
%
% A tomArray is a multidimensional array that contains data either in the
% shape of a numeric arrray, or a tomsym object.
%
% TomArrays makes it possible to work with multidimensional tomsym
% expressions in the same way as ordinary Matlab arrays work.
% A number of Matlab operators and functions are overloaded to work with
% tomArrays, including +, -, .*, ./, sum, permute, reshape, squeeze, ...
%
% In addition to working like ordinary multidimensional matlab arrays,
% tomArrays can also be used with "named indexes". This is achieved by
% using a string constant as an index into the array. When operators such
% as + or .* are used on tomArrays with named indexes, the arrays will
% automatically be permuted so that the indexes with the same name line up.
% If there are indexes that are used by only one of the operand, then a
% repmat operation will be applied to the other in order to make both
% arrays the same size.
%
% For example, assuming that A is a 2-by-3 tomArray and B is a 3-by-2
% tomArray, it is possible to do the following:
% 
%   i='i'; j='j'; k='k'; % Shorthand notation for string constants.
%
%   A(i,1) % Returns a 1-D tomArray of length 3 and index name 'i'.
%   A(:,1) % The same tomArray as above, but without index name.
%   A(i,j)+B(j,i)  % A 2-by-3 tomArray with index names 'i', 'j'.
%   A(i,j).*B(j,k) % A 2-by-3-by-2 tomArray indexed 'i', 'j', 'k'.
%   sum(A(i,j).*B(j,k),j) % A 2-by-2 tomArray (the matrix product).
%
% Making a singleton slice in a tomArray (as in A('i',1)) reduces the
% number of dimensions, and of course the removed index doesn't need a
% name. However, when making a wider slice, it is necessary to either name
% that index as well, or leave the entire array unnamed. For example:
%
%    % Same shorthand notation as above: i='i', etc.
%    A(i,{j 1:end-1}) % Returns a 2-by-2 tomArray indexed 'i', 'j'
%    A(:, 1:end-1)  % The same tomArray as above, but without index names.
%
% Note that, since i, j, etc. are character constants, AMPL-style indexing
% A(i,j+1) does not work as expected. (Instead j+1 will evaluate to 107,
% since the ASCII code for 'j' is 106.) In this case it is better to use
% the tomArrayIdx class for indexing. Example:
%
%    j = tomArrayIdx('j', 1:2);
%    A(i,j+1) % Returns a 2-by-2 tomArray indexed 'i', 'j'
%
% A tomArray can be converted back into its underlying data format
% (typically a tomSym expression or a numeric array) using the command
% unArray.
%
% It is also possible to convert any tomArray back into a tomSym or numeric
% vectors using the vec() function, or using the syntax: A(:)
%
% See also: tomArrayIdx

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-26 by rutquist for TOMLAB release 7.7

if nargin==0
    M = [];
end

if nargin==1
    if isa(X,'tomArray')
        M = X;
    elseif isnumeric(X) || isa(X,'tomSym')
        M = struct;
        M.X = X;
        if length(size(X))==2 && size(X,2)==1
            M.sz = size(X,1);
        elseif length(size(X))==2 && size(X,1)==1
            M.sz = size(X,2);
        else
            M.sz = size(X);
        end
        M.ni = {};
        M = class(M,'tomArray');
    elseif isa(X,'tomArrayIdx');
        M = struct;
        M.X = double(X);
        M.sz = length(M.X);
        M.ni = {char(X)};
        M = class(M,'tomArray');
    else
        error(['Class ' class(X) 'cannot be converted into a tomArray.']);
    end
else
    if ischar(X)
        if nargin>=4
            X = tom(X,prod(sz),1,intflag);
        else
            X = tom(X,prod(sz),1);
        end
    end
    
    sz = sz(:)';
    if ~numel(X) == prod(sz)
        error('Number of elements must match specified size');
    end

    if isa(X,'tomArray')
        X = X.X;
    end
    
    M = struct;
    M.X = X;
    M.sz = sz;
    if nargin<3 || isempty(ni);
        M.ni = {};
    else
        if ~iscell(ni)
            error('Index list must be a cell arrary.');
        end
        M.ni = ni;
    end
    M = class(M,'tomArray');
    checkIndexes(M);
end

% Establish priority over tomSym
superiorto('tomSym');
