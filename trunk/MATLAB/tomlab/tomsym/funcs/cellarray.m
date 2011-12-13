function y = cellarray(varargin)
% Y = CELLARRRAY(ITEM, ....) put arguments into a cell array.
%
% This function is used by tomSym to handle symbolic arrays. It should
% never need to be called explicitly.
%
% Note: symbolic cell arrays are only useable under special circumstances.

y = {varargin{:}};
