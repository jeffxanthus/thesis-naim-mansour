% function F0 = multiMin_1_par(X0,Prob)
%
% Compute loop 1 in multiMin with parfor using Prob.Threads threads

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written July 13, 2006.    Last modified July 23, 2011.

function F0 = multiMin_1_par(X0,Prob)
M  = size(X0,2);
F0 = zeros(1,M);
for P = 1:M
    F0(1,P) = nlp_f(X0(:,P), Prob);         % Compute f(x_0) for each x_0
end

% MODIFICATION LOG:
%
% 110723  hkh  Written
