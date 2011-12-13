% lp2tp Convert Transportation LP to transportation formulation
%       defining supply, demand and Cost (distance) matrix
%
% function [s, d, C] = lp2tp(A,b,c)
%
%                    m = M + N; n = M * N
%
%   Input:
%           A 	     Node-arc incidence matrix. A is m x n.
%           b 	     Supply/demand vector. Length m.
%           c 	     Cost vector. Length n.
%
%   Output:
%           s        Supply, M by 1
%           d        Demand, N by 1
%           C        Cost matrix, M by N

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomopt.com
% Copyright (c) 2003-2004 by Tomlab Optimization Inc., $Release: 4.0.6$
% Written Feb 26, 2003.   Last modified Mar 1, 2003.

function [s, d, C] = lp2tp(A,b,c)

i = find(A(1,:) ==0);
n = i(1) - 1;
if n <= 1
    i = find(A(1,:) == 1);
    n = median(i(2:end)-i(1:end-1));
    m = length(b) - n;
    b = b(:);
    d = b(1:n);
    s = b(n+1:n+m);
else
    m = length(b) - n;
    b = b(:);
    s = b(1:m);
    d = b(m+1:m+n);
end

C = reshape(c(:),n,m)';

% MODIFICATION LOG
%
% 030226 hkh Written
% 030301 hkh Improve help, use C for cost matrix, transpose of C