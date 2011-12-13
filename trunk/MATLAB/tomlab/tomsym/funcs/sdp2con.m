function [c,x0] = sdp2con(m,x0,sname,issym) 
% sdp2con - Convert a semidefinite constraint into a noninear one
%
% c = sdp2con(m) creates a set of nonlinear constraints that can only be
% satisfied if m is a positive definite matrix. This makes it possible to
% use constraints from semidefinite programming with ordinary NLP solvers.
%
% The input m must be a symmetric tomSym matrix.
%
% [c,x0] = sdp2con(m,x0,sname,issym) specifies the following optional inputs:
%
%   - x0 is a struct representing an initial guess. An initial guess for the
%     new unknowns will be appended to x0.
%   - sname specifies the name of the new symbol that will be introduced
%   - issym if true indicates that the matrix m is symmetric.
%
% See also: MI

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2011 by Tomlab Optimization Inc.
% $Id$

if nargin < 2
    x0 = struct;
end

if nargin < 3
    sname = '';
end

if tomCmp(m,'positiveSemidefinite')
    m = operand(1,m);
end

if nargin < 4
    issym = isequal(m,m');
end

if ~isa(m,'tomSym') || size(m,1) ~= size(m,2)
    error('input must be a symbolic square matrix');
end

P = pattern(m);
if issym && ~isequal(P,P')
    error('MI sparsity pattern does not appear to be symmetric.');
end

N = size(P,1);

s = symamd(P,[-1,0]);
A = P(s,s);
if ~issym
    A = A | A';
end

L = false(N,N); 
for j=1:N
    L(j,j) = A(j,j) | any(L(j,1:j-1)); 
    for i=j+1:N
        L(i,j) = A(i,j) | any(L(i,1:j-1) & L(j,1:j-1));
    end
end

[ix,jx] = find(L);

symb = tom(sname,length(ix),1);

x = sparse(ix,jx,symb,N,N);

if issym
    ms = m(s,s);
else
    ms =  0.5*(m(s,s)+m(s,s)');
end

ixr = tril(double(L)*double(L)' ~= 0);

c = {ms(ixr) == lookup(x*x',ixr), symb(ix==jx) >= 0 };

if nargout >= 2
    if ~isstruct(x0)
        x0 = tom2struct(x0);
    end
    M = subs(ms,x0);
    if isnumeric(M)
        if issym && any(abs(vec(M-M'))>=4*eps(max(abs(M(:)))))
            error('MI does not appear to be symmetric when guess is applied');
        end
        l = eigs(M,1,'sa');
        if l<0
            M = M-1.01*l*eye(N);
        end
        L0 = chol(M,'lower');
        x0.(char(symb)) = L0(L(:));
    else
        x0.(char(symb)) = repmat(1/length(symb),size(symb));
    end
end
