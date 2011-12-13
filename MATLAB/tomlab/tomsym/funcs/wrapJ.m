function J = wrapJ(info,m,varargin)
% wrapJ - Compute the Jacobian of a wrapped function.
%
% J = wrapJ(info,M,...) returns the Jacobian of info.fun with respect to
% the M:th input argument
%
% J will be a Jacobian matrix on the form that is used by tomSym.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-10-02 by rutquist for TOMLAB release 7.7

if ~isfield(info,'nOut')
    info.nOut = 1;
end

if isfield(info,'Jfuns') && length(info.Jfuns)>=m;
    jfun = info.Jfuns{m};
else
    jfun = 'MAD';
end

if ~isfield(info,'id')
    info.id = char('A'+floor(25*rand(1,14)));
end

switch lower(jfun)
    case 'mad'
        J = madWrap(info.nOut,info.fun,m,varargin{:});
    case 'fdjac'
        J = FDJacR(info.nOut,info.fun,m,...
            info.Jpatterns{m},info.Jcolidx{m},[info.id num2str(m)],...
            varargin{:});
    case 'cdjac'
        if isfield(info,'CDh')
            CDh = info.CDh{m};
        else
            CDh = 0;
        end
        J = CDJac(info.nOut,info.fun,m,...
            CDh,info.Jpatterns{m},info.Jcolidx{m},...
            varargin{:});
    otherwise
        J = feval(jfun,varargin{:});
end

