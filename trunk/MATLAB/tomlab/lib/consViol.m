% function [h_L1, h] = consViol(z,L,U,absviol)
%
% Compute constraint violation

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written July 13, 2006.    Last modified July 23, 2011.

% -----------------------------------------------
function [h_L1, h] = consViol(z,L,U,absviol)
% -----------------------------------------------
if absviol == 1
   h =-min(0,z-L)-min(0,U-z);
else
   h =-min(0,(z-L)./max(1,abs(L)))-min(0,(U-z)./max(1,abs(U)));
end
h_L1 = sum(h);

%fprintf('h_L1  %25.10f ',h_L1);
%h_L2 = norm(h);
%fprintf('h_L2  %25.10f\n',h_L2);


% MODIFICATION LOG:
%
% 110723  hkh  Written
