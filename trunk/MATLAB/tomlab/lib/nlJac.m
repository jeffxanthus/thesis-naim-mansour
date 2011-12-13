% nlJac.m
%
% function [J, dg]=nlJac(x, Prob)
%
% nlJac calls the TOMLAB gateway functions nlp_g and nlp_dc,
% to evaluate the Jacobian of the objective and constraint function.
%
% The output is divided into inequality and equality constraints
% and transformed to fit the formulation:
%
%  g(x) == 0  (linear constraints first, then nonlinear)
%  g(x) >= 0  (linear constraints first, then nonlinear)
%
% nlJac is used as callback from Schittkowski routines

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2003 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written May 14, 2003.   Last modified May 14, 2003.

function [J, dg]=nlJac(x, Prob)

J  = full(nlp_J( x, Prob) );
dc = nlp_dc( x, Prob);

% Allocate full constraint Jacobian matrix

m = Prob.m;
dg = zeros(m,Prob.N);

mEQ = Prob.mEQ;
mIN = Prob.mIN;

AixEQ  = Prob.AixEQ;
AixLow = Prob.AixLow;
AixUpp = Prob.AixUpp;

cixEQ  = Prob.cixEQ;
cixLow = Prob.cixLow;
cixUpp = Prob.cixUpp;

% First compute Jacobian of equality constraints

if mEQ > 0
   if isempty(AixEQ)
      dg(1:mEQ,:) = dc(cixEQ,:);
   elseif isempty(cixEQ)
      dg(1:mEQ,:) = Prob.A(AixEQ,:);
   else
      dg(1:length(AixEQ),:)     = Prob.A(AixEQ,:);
      dg(length(AixEQ)+1:mEQ,:) = dc(cixEQ,:);
   end
end
% Compute Jacobian of inequality constraints
if mIN > 0
   if ~(isempty(AixLow) && isempty(AixUpp))
      if isempty(AixLow)
         m2 = mEQ+length(AixUpp);
         dg(mEQ+1:m2,:) =  -Prob.A(AixUpp,:);
      elseif isempty(AixUpp)
         m2 = mEQ+length(AixLow);
         dg(mEQ+1:m2,:) = Prob.A(AixLow,:);
      else
         m2 = mEQ+length(AixLow);
         dg(mEQ+1:m2,:) = Prob.A(AixLow,:);
         dg(m2+1:m2 + length(AixUpp),:) = -Prob.A(AixUpp,:);
         m2 = m2+length(AixUpp);
      end
   else
      m2 = mEQ;
   end
   if m2 ~= m
      if isempty(cixLow)
         m3 = m2+length(cixUpp);
         dg(m2+1:m3,:) = -dc(cixUpp,:);
      elseif isempty(AixUpp)
         m3 = m2+length(cixLow);
         dg(m2+1:m3,:) = dc(cixLow,:);
      else
         m3 = m2+length(cixLow);
         dg(m2+1:m3,:) = dc(cixLow,:);
         dg(m3+1:m,:) = -dc(cixUpp,:);
      end
   end
end

% MODIFICATION LOG:
%
% 030514  hkh  Written, based on nlgrad