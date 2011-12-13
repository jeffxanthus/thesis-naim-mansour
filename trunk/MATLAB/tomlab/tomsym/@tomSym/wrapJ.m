function y = wrapJ(info,m,varargin)
% tomSym/wrapJ - Overloaded

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if ~isfield(info,'id')
    info.id = char('A'+floor(25*rand(1,14)));
end

y = tomSym(mfilename,info.sz1*info.sz2,numel(varargin{m}),info,m,varargin{:});
