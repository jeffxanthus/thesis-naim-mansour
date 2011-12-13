% function Prob = schedulingofteleviasatellite(traffic)
%
% Creates a TQBS for scheduling of telecommunications via satellite
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Nov 22, 2005.   Last modified Nov 22, 2005.

function TQBS = schedulingofteleviasatellite(traffic)

if nargin < 1
   error('The function requires 1 inputs');
end

if isempty(traffic)
   error('One of the inputs are empty');
end

% Convert the traffic matrix so all rows and columns equal maxsum
n1 = size(traffic, 1);
rowsum = sum(traffic,1);
colsum = sum(traffic,2);
LB = max(max([rowsum;colsum']));
q = zeros(n1,n1);
for t=1:n1
   for r=1:n1
      q = min([LB-rowsum(t);LB-colsum(r)]);
      traffic(r,t) = traffic(r,t) + q;
      rowsum(t) = rowsum(t) + q;
      colsum(r) = colsum(r) + q;
   end
end

TQBS = traffic;

% MODIFICATION LOG
%
% 051122 med   Created.