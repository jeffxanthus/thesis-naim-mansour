% function [mA, Ax, bEqual, b_L, b_U, A] = LinearConstr(Prob)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written June 28, 1999. Last modified Feb 17, 2000.

function [mA, Ax, bEqual, b_L, b_U, A] = LinearConstr(Prob)

if nargin < 1
   error('LinearConstr must have 1 parameter Prob as input');
end

% Linear constraints

b_L=  Prob.b_L(:);
b_U=  Prob.b_U(:);
A=    Prob.A;

if isempty(A)
   mA=0;
else
   mA=size(A,1);
end

if mA > 0
   if isempty(b_U),b_U= Inf*ones(mA,1); end
   if isempty(b_L),b_L=-Inf*ones(mA,1); end

   if isempty(Prob.x_0)
      Ax=[];
   else
      Ax=A*Prob.x_0;
   end
   bEqual=eq(b_L,b_U);
else 
   Ax=zeros(0,1); % Because b_L(:) makes empty vector be size (0,1)
   bEqual=zeros(0,1);
end