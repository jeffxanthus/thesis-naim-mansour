function y = norm(a,P)
% tomSym/norm - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-10-28 by rutquist for TOMLAB release 7.7

if nargin<2
    P=2;
end

if size(a,1)==1 || size(a,2)==1
    % Norm of a vector
    switch P
        case 1
            y = sum(abs(a));
        case 2
            y = sqrt(sum(a.^2));
        case inf
            y = max(abs(a));
        case -inf
            y = min(abs(a));
        otherwise
            y = sum(abs(a).^P)^(1/P);
    end
else
    % Norm of a matrix
    switch P
        case 1
            y = max(sum(abs(a)));
        case 2
            info = struct;
            info.fun = 'norm';
            info.n = 1;
            info.sz1 = 1;
            info.sz2 = 1;
            info.Jpatterns = {ones(1,numel(a))};
            info.Jfuns = {'CDJac'};
            info.CDh = {1e-6};
            y = wrap(info,a,2);
        case inf
            y = max(sum(abs(a),2));
        case 'fro'
            y = sqrt(sum(diag(a'*a)));
    end
end
