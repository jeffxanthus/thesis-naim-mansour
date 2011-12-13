% ksDef.m
%
% function Prob=ksDef(Prob,sense)
%
% ksDef creates index vectors and counters to be used in the
% Klaus Schittkowski (ks) routines.
%
% The structure information is passed to routines nlfunc and nlgrad.
%
% The Tomlab constraints
% b_L <= A * x <= b_U
% c_L <= c(x)  <= c_U
%
% are transformed to the ks format.
%
% The output is divided into inequality and equality constraints
% and transformed to fit the formulation:
%
%  g(x) == 0  (linear constraints first, then nonlinear)
%  g(x) >= 0  (linear constraints first, then nonlinear)
%
% After conversion the constraints are
%
%  A * x - b_L    = 0
%  c(x)  - c_L    = 0
%
%  A * x - b_L   >= 0
%  b_U   - A * x >= 0
%  c(x)  - c_L   >= 0
%  c_U   - c(x)  >= 0
%
% There are Prob.mEQ equality   constraints
% There are Prob.mIN inequality constraints
% The total number of constraints after conversion is Prob.m
%
% The fields set to be used in Prob are:
% AixEQ, AixLow, AixUpp
% cixEQ, cixLow, cixUpp
% m, mEQ, mIN
%
% The 'sense' parameter is used to flip the inequalities
% from '>=' to '<='. Set to -1 for '<='. Set to 1,empty or omit the
% 'sense' parameter entirely for default '>='.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2003-2009 by Tomlab Optimization Inc., Sweden. $Release: 7.3.0$
% Written April 27, 2003.   Last modified Aug 13, 2009.

function Prob = ksDef(Prob,sense)

if nargin < 2
   % Default KS >= sense for single-sided constraints
   sense = 1;
end

if ~isempty(Prob.A)
   if isempty(Prob.b_L)
      Prob.AixUpp = find(~isinf(Prob.b_U));
      Prob.AixLow = [];
      Prob.AixEQ  = [];
   elseif isempty(Prob.b_U)
      Prob.AixLow = find(~isinf(Prob.b_L));
      Prob.AixUpp = [];
      Prob.AixEQ  = [];
   else
      ixEQ        = Prob.b_L==Prob.b_U & ~isinf(Prob.b_L);
      Prob.AixLow = find(~isinf(Prob.b_L) & ~ixEQ);
      Prob.AixUpp = find(~isinf(Prob.b_U) & ~ixEQ);
      Prob.AixEQ  = find(ixEQ);
   end
else
   Prob.AixEQ  = [];
   Prob.AixLow = [];
   Prob.AixUpp = [];
end

if ~(isempty(Prob.c_L) && isempty(Prob.c_U))
   if isempty(Prob.c_L)
      Prob.cixUpp = find(~isinf(Prob.c_U));
      Prob.cixLow = [];
      Prob.cixEQ  = [];
   elseif isempty(Prob.c_U)
      Prob.cixLow = find(~isinf(Prob.c_L));
      Prob.cixUpp = [];
      Prob.cixEQ  = [];
   else
      ixEQ        = Prob.c_L==Prob.c_U & ~isinf(Prob.c_L);
      Prob.cixLow = find(~isinf(Prob.c_L) & ~ixEQ);
      Prob.cixUpp = find(~isinf(Prob.c_U) & ~ixEQ);
      Prob.cixEQ  = find(ixEQ);
   end
else
   Prob.cixEQ  = [];
   Prob.cixLow = [];
   Prob.cixUpp = [];
end

Prob.mEQ = length(Prob.AixEQ) + length(Prob.cixEQ);
Prob.mIN = length(Prob.AixLow) + length(Prob.AixUpp) + ...
           length(Prob.cixLow) + length(Prob.cixUpp);
Prob.m   = Prob.mEQ + Prob.mIN;
Prob.sense = sense;

% MODIFICATION LOG:
%
% 030427  hkh  Written
% 040109  ango Expanded functionality with 'sense' input
% 090813  med  mlint check