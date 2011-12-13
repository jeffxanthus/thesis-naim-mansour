function y = ppnval(PPN, varargin)
% TOMSYM/PPNVAL - Overloaded function

y = tomSym(mfilename, size(varargin{1},1), size(varargin{1},2), PPN, varargin{:});
