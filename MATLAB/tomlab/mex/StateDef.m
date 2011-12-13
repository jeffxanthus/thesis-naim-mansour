% StateDef defines the state variables in the Result structure
%
% function Result = StateDef(Result,x_k,Ax,c_k,xTol,bTol,cTol,bl,bu,Order)
%
% INPUT:
% Result The current Result structure, where fields are to added
% x_k    The variables at the optimal solution
% Ax     The linear constraints evaluated at x_k, i.e. A * x_k
% c_k    The nonlinear constraints evaluated at x_k
% xTol   Variable tolerance
% bTol   Linear constraint tolerance
% cTol   Nonlinear constraint tolerance
% bl     Lower bounds on variables, nonlinear and linear constraints
% bu     Upper bounds on variables, nonlinear and linear constraints
% Order  If > 0 the order of bl,bu is assumed to be
%        (variables, linear constraints, nonlinear constraints)
%        Default: Order = 0
%
% OUTPUT:
% Result   Structure with results (see ResultDef.m):
%          Fields xState, bState, cState and QP.B set

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written July 24, 2000. Last modified Jun 6, 2008.

function Result = StateDef(Result,x_k,Ax,c_k,xTol,bTol,cTol,bl,bu,Order)

if nargin < 10, Order = []; end
if isempty(Order), Order = 0; end

% State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
n  = length(x_k);
m1 = length(Ax);
m2 = length(c_k);

Result.xState=double(x_k-xTol <= bl(1:n)) + 2*double(x_k >= bu(1:n)-xTol);

if m1 > 0
   if Order > 0
      k=n;
   else
      k=n+m2;
   end
   Result.bState=double(Ax-bTol<=bl(k+1:k+m1))+2*double(Ax>=bu(k+1:k+m1)-bTol);
end
if ~isempty(c_k)
   if Order > 0
      k=n+m1;
   else
      k=n;
   end
   Result.cState=double(c_k-cTol<=bl(k+1:k+m2))+2*double(c_k>=bu(k+1:k+m2)-cTol);
end

B = zeros(n,1);

B(find( ~( (x_k > bl(1:n) + xTol) & (x_k < bu(1:n) - xTol) & ...
   (bl(1:n)~=bu(1:n))))) = 1;
B(find(x_k >= bu(1:n) - xTol)) = -1;
Result.QP.B=B;

% MODIFICATION LOG:
%
% 000724 hkh  Written
% 020819 hkh  Revised for Matlab 6.5 logical handling
% 021022 hkh  Use double(zeros( )) to make it work for Matlab 6.5
% 030902 ango Went back to (revised version of) old code for calculating
%             xState,bState,cState
% 041201 hkh  Add input Order and allow different order of constraints
% 041201 hkh  Correct and revise help
% 080606 med  Cleaned