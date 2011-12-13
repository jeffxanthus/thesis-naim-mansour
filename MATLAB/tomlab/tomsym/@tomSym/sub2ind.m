function ind = sub2ind(siz,i,j,varargin)
% tomSym/sub2ind - Overloaded function.
%
% sub2ind(siz,i,j) returns siz(1)*(j-1) + i
%
% This is a simplified version of Matlab's sub2ind, which only works for
% two-dimensional indexes.
%
% It is provided so that functions that rely on sub2ind will work with
% tomSym symbolic indexes.

if length(siz)~=2 || ~isempty(varargin)
    error('sub2ind with tomSym arguments only works with 2D arrays.');
end

ind = lookup(siz,1)*(j-1) + i;
