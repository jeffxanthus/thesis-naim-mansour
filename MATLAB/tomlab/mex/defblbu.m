% Define lower and upper bound arrays for the SOL family of solvers
%
% Inf are changed to BIG (=1E20), -Inf to -BIG.
%
% function [bl, bu, n, m1, m2] = defblbu(Prob, BIG, Order);
%
% INPUT:
%
%   Prob     Used fields in structure Prob:
%            x_L    Lower bounds on x
%            x_U    Upper bounds on x
%            b_L    Lower bounds on linear constraints
%            b_U    Upper bounds on linear constraints
%            c_L    Lower bounds on nonlinear constraints
%            c_U    Upper bounds on nonlinear constraints
%
%   BIG      Numeric Inf value. Default 1E20
%   Order    The order of the constraints in output bl and bu
%            If Order > 0 (default)
%            (x, linear constraints, nonlinear constraints)
%            otherwise
%            (x, nonlinear constraints, linear constraints)
%
% OUTPUT:
%   bl       Lower bound vector, input to SOL solvers
%   bu       Upper bound vector, input to SOL solvers
%            The order is (x, linear constraints, nonlinear constraints)
%   n        Number of variables
%   m1       Number of linear constraints
%   m2       Number of nonlinear constraints

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Mar 2, 1998.    Last modified May 26, 2007.

function [bl, bu, n, m1, m2] = defblbu(Prob, BIG, Order)

if nargin < 3
   Order = [];
   if nargin < 2
      BIG=1E20;
end, end

if isempty(Order), Order = 1; end

% Simple bounds
n = Prob.N;

% Linear constraints
% m1=max([length(Prob.b_L),length(Prob.b_U),size(Prob.A,1)]);
m1 = Prob.mLin;

% Nonlinear constraints
% m2=max(length(Prob.c_L),length(Prob.c_U));
m2 = Prob.mNonLin;

%n=max([length(bl),length(bu),length(Prob.x_0)]);
%if n==0, error('defblbu: cannot determine n - b_L, b_U and x_O are all empty');end
%

bl=zeros(n+m1+m2,1);
bu=zeros(n+m1+m2,1);

nL = length(Prob.x_L);
nU = length(Prob.x_U);

if nL >= n        
   bl(1:n)    = Prob.x_L(1:n);
elseif nL == 0        
   % Nothing given - set all -BIG
   bl(1:n)    = -BIG;
else
   % Pad missing elements with -BIG
   warning('defblbu: padding x_L with -BIG')
   bl(1:nL)   = Prob.x_L;
   bl(nL+1:n) = -BIG;
end
bl=bl(:);

if nU >= n        
   bu(1:n)    = Prob.x_U(1:n);
elseif nU == 0        
   % Nothing given - set all BIG
   bu(1:n)    = BIG;
else
   % Pad missing elements with -BIG
   warning('defblbu: padding x_U with BIG')
   bu(1:nU)   = Prob.x_U;
   bu(nU+1:n) = BIG;
end
bu=bu(:);

if m1 > 0
   if Order > 0
      k = n;
   else
      k = n+m2;
   end
   m1L = length(Prob.b_L);
   m1U = length(Prob.b_U);

   if m1L >= m1
      bl(k+1:k+m1)     = Prob.b_L;
   elseif m1L == 0
      bl(k+1:k+m1)     = -BIG;
   else
      warning('defblbu: padding b_L with -BIG')
      bl(k+1:k+m1L)    = Prob.b_L;
      bl(k+m1L+1:k+m1) = -BIG;
   end
   if m1U >= m1
      bu(k+1:k+m1)     = Prob.b_U;
   elseif m1U == 0
      bu(k+1:k+m1)     = BIG;
   else
      warning('defblbu: padding b_U with BIG')
      bu(k+1:k+m1U)    = Prob.b_U;
      bu(k+m1U+1:k+m1) = BIG;
   end
end

if m2 > 0
   if Order > 0
      k = n+m1;
   else
      k = n;
   end
   m2L = length(Prob.c_L);
   m2U = length(Prob.c_U);
   if m2L >= m2
      bl(k+1:k+m2)     = Prob.c_L;
   elseif m2L == 0
      bl(k+1:k+m2)     = -BIG;
   else
      warning('defblbu: padding c_L with -BIG')
      bl(k+1:k+m2L)    = Prob.c_L;
      bl(k+m2L+1:k+m2) = -BIG;
   end
   if m2U >= m2
      bu(k+1:k+m2)     = Prob.c_U;
   elseif m2U == 0
      bu(k+1:k+m2)     = BIG;
   else
      warning('defblbu: padding c_U with BIG')
      bu(k+1:k+m2U)    = Prob.c_U;
      bu(k+m2U+1:k+m2) = BIG;
   end
end

if ~isinf(BIG)
   % Adjust everything to be in the interval [-BIG, BIG]
   ix=find(bl < -BIG);
   bl(ix)=-BIG;

   ix=find(bu > BIG);
   bu(ix)=BIG;
end

% MODIFICATION LOG
%
%  980921  hkh  Changed prob to Prob. Changed comments
%  990916  hkh  Compute m1 and m2. More safe computation for Mideva.
%  020920  hkh  Using m1 and not m2 in second part
%  030129  ango Pad with (+/-)BIG if too few (but >0) values given
%  031114  hkh  Change default BIG from 1E20 to 1E12
%  040109  hkh  Use mLin, mNonLin instead of computing m1,m2
%  041202  hkh  Modify and correct help, make order of bl/bu clear
%  041202  hkh  Add parameter Order, if > 0, reverse order
%  041202  hkh  Total revision, avoid vector expansion
%  041202  hkh  Only test using BIG if ~isinf