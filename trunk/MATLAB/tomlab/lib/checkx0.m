% Check if x_0 is inside bounds. If not call pbuild to initiate
% savings of search directions.
%
% function checkx0(x_0, x_L, x_U);
%
% This routine is used for routines which first adjust the starting value
% into the feasible region before calling the function value routine
%
% x_0      Starting values
% x_L      Lower bounds on x
% x_U      Upper bounds on x

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Nov 20, 1999.  Last modified June 22, 1999.

function checkx0(x_0, x_L, x_U)

if isempty(x_0)
   return
end

x_s=x_0;

if ~isempty(x_U)
   % Fix that x is below upper bounds
   x_s(1:length(x_U))=min(x_0(1:length(x_U)),x_U); 
end

if ~isempty(x_L)
   % Fix that x is greater or equal to lower bounds
   x_s(1:length(x_L))=max(x_s(1:length(x_L)),x_L); 
end

if sum(abs(x_s-x_0)) > 0 % Moved x_0 inside bounds
   pbuild(x_0,NaN);
end

% MODIFICATION LOG:
%
% 981107  hkh  Check if x_0 is empty, then return