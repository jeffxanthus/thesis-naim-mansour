% ====================================================================
% function onB = nOnBound(X,x_L,x_U,epsX)
% ====================================================================
% Find number of components of X that is on or very close to the simple bounds
%
% epsX is a tolerance given as input, default 1.7E-6
%
% onB is computed as:
%
% sum(abs(x_L-X)<=epsX*max(1,abs(x_L)) | abs(X_U-x)<=epsX*max(1,abs(x_U)));
%
% If X is matrix, a vector of onB values are returned

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2006 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Nov 24, 2006. Last modified Oct 20, 2009.

function onB = nOnBound(X,x_L,x_U,epsX)
if nargin < 4
   epsX = 1.7E-6;
end
m = size(X,2);
if m == 1
   onB = sum(abs(x_L-X)<=epsX*max(1,abs(x_L)) | abs(x_U-X)<=epsX*max(1,abs(x_U)));
else
   onB = zeros(1,m);
   for i=1:m
       onB(i) = sum(abs(x_L-X(:,i))<=epsX*max(1,abs(x_L)) | ...
                    abs(x_U-X(:,i))<=epsX*max(1,abs(x_U)));
   end
end
% xprinti(onB,'nOnBound:')

% MODIFICATION LOG
%
% 061124  hkh  First version written
% 091020  hkh  If X matrix, return vector of onB values
