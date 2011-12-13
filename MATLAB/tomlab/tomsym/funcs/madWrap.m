function J = madWrap(m,fun,n,varargin)
% madWrap - Compute a Jacobian using MAD.
%
% J = madWrap(FUN,N,...) uses MAD to call FUN(...) and returns the Jacobian
% matrix with respect to the N:th input argument.
%
% FUN must be the name of an existing function.
% N must be an integer between one and nargin(FUN).
%
% J = madWrap(M,FUN,N,...) computes the Jacobian matrix of the Mth output
% argument, instead of the first one.
%
% J will be a Jacobian matrix on the form that is used by tomSym.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if ischar(m)
    varargin = {n, varargin{:}};
    n = fun;
    fun = m;
    m = 1;
end

% Setup input arguments
args = varargin;
an = varargin{n};
nd = numel(an);
da = zeros([size(an) nd]);
da((nd+1)*(0:nd-1)+1) = 1;
args{n} = fmad(an,da);

% Evaluate function and derivative
y = cell(1,m);
[y{:}] = feval(fun,args{:});
y = y{m};

% Reshape MAD derivatives into tomSym Jacobian.
J = reshape(getderivs(y),[getvalue(numel(y)),nd]);
