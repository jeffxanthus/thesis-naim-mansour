% function beta = convrate(m)
%
% Compute estimate of convergence rate in each step from global matrix p_dx
%
% ||x_k+1 - x_k   ||       ||p_k+1 ||
% ------------------   =   ----------    = beta_k
% ||x_k   - x_k-1 ||^m     ||p_k   ||^m
%
% where m is the order of convergence. default m=1.
%
% beta_1 == NaN
%
% Globals:
%   p_dx:	Matrix with all the search directions
%   alphaV: Vector with the alfa steps (line search lengths)
%

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1997-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Nov 11, 1997. Last modified June 22, 1999.

function beta = convrate(m)

if nargin < 1, m=[]; end
if isempty(m), m=1; end
if m <= 0, m=1; end

beta=NaN;

global p_dx

if isempty(p_dx), return, end

n=size(p_dx,2);

if n==0, return, end

beta=NaN*ones(n,1);

pk=norm(p_dx(:,1));

for i=2:n
    pk1=norm(p_dx(:,i));
    if pk~=0
        if m==1
            beta(i)=pk1/pk;
        else
            beta(i)=pk1/pk^m;
        end
    end
    pk=pk1;
end