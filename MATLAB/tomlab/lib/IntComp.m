% -------------------------------------------------
% Check if the x, given as columns in X, is integer valued
%
% -------------------------------------------------
% function [IC,Iidx] = IntComp(X,IntVars,eps_I)
% -------------------------------------------------
%
% INPUT:
%   X       Columns of x vectors
%   IntVars Indices for the integer variables
%   eps_I   Tolerance how close to integer is OK
% OUTPUT:
%   IC      Number of integer components in each column in X
%   IIdx    For each component in X, 1 = integer valued, 0 = not integer

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Oct 15, 2009.   Last modified Oct 18, 2009.

function [IC,Iidx] = IntComp(X,IntVars,eps_I)

m  = size(X,2);
IC = zeros(1,m);
if isempty(IntVars)
   Iidx = zeros(0,m);
else
   Iidx = zeros(length(IntVars),m);
   for i=1:m
       z         = round(X(IntVars,i));
       Iidx(:,i) = abs(z-X(IntVars,i)) <= eps_I;
       IC(i)     = sum(Iidx(:,i));
   end
end
Iidx = logical(Iidx);

% MODIFICATION LOG
%
% 091015  hkh  Written
% 091018  hkh  Set Iidx logical
