% function PrintTSP(Result)
%
% Prints TSP (Travelling Salesman Problem) results.
%
% INPUT PARAMETERS
% Result       TOMLAB Result structure

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 5.4.0$
% Written Jul 31, 2006.   Last modified Jul 31, 2006.

function PrintTSP(Result)

if nargin < 1
    error('PrintTSP requires the Result as input');
end

a     =  size(Result.Prob.TSP.C,1);
temp  = reshape(Result.x_k(1:a^2),a,a);
route = zeros(a+1,1);
route(1,1) = 1;
for i=1:a
    route(i+1,1) = find(temp(route(i),:) > 0.5);
end

disp(['Shortest route is: ' num2str(route')]);

% MODIFICATION LOG
%
% 060731  med  Created