% nlfunc.m
%
% function [f, g]=nlfunc(x, Prob)
%
% nlfunc calls the TOMLAB gateway functions nlp_f and nlp_c,
% to evaluate the objective and constraint function.
%
% The output is divided into inequality and equality constraints
% and transformed to fit the formulation:
%
%  g(x) == 0  (linear constraints first, then nonlinear)
%  g(x) >= 0  (linear constraints first, then nonlinear)
%
% nlfunc is used as callback from Schittkowski routines and GENO. 

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2003-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.5.0$
% Written April 27, 2003.   Last modified Sep 4, 2006.

function [f, g] = nlfunc(x, Prob)

f = nlp_f( x, Prob);

% Allocate full constraint vector
m = Prob.m;
g = zeros(m,1);

% Exit early if completely unconstrained
if m==0
   return
end

% If only linear constraints, do not call nlp_c
if Prob.mNonLin>0
   c = nlp_c( x, Prob);
else
   c = [];
end

% First compute equality constraints
mEQ = Prob.mEQ;
mIN = Prob.mIN;

AixEQ  = Prob.AixEQ;
AixLow = Prob.AixLow;
AixUpp = Prob.AixUpp;

cixEQ  = Prob.cixEQ;
cixLow = Prob.cixLow;
cixUpp = Prob.cixUpp;

if mEQ > 0
   if isempty(AixEQ)
      g(1:mEQ) = c(cixEQ) - Prob.c_L(cixEQ) ;
   elseif isempty(cixEQ)
      g(1:mEQ) = Prob.A(AixEQ,:) * x - Prob.b_L(AixEQ);
   else
      g(1:length(AixEQ))     = Prob.A(AixEQ,:) * x - Prob.b_L(AixEQ);
      g(length(AixEQ)+1:mEQ) = c(cixEQ) - Prob.c_L(cixEQ) ;
   end
end
% Compute inequality constraints
sense=Prob.sense(1); 
if mIN > 0
   if ~(isempty(AixLow) && isempty(AixUpp))
      if isempty(AixLow)
         m2 = mEQ+length(AixUpp);
         g(mEQ+1:m2) = sense*(Prob.b_U(AixUpp) - Prob.A(AixUpp,:) * x);
      elseif isempty(AixUpp)
         m2 = mEQ+length(AixLow);
         g(mEQ+1:m2) = sense*(Prob.A(AixLow,:) * x - Prob.b_L(AixLow));
      else
         m2 = mEQ+length(AixLow);
         g(mEQ+1:m2) = sense*(Prob.A(AixLow,:) * x - Prob.b_L(AixLow));
         g(m2+1:m2 + length(AixUpp)) = sense*(Prob.b_U(AixUpp) - Prob.A(AixUpp,:) * x);
         m2 = m2+length(AixUpp);
      end
   else
      m2 = mEQ;
   end
   if m2 ~= m
      if isempty(cixLow)
         m3 = m2+length(cixUpp);
         g(m2+1:m3) = sense*(Prob.c_U(cixUpp) - c(cixUpp));
      elseif isempty(cixUpp)
         m3 = m2+length(cixLow);
         g(m2+1:m3) = sense*(c(cixLow) - Prob.c_L(cixLow));
      else
         m3 = m2+length(cixLow);
         g(m2+1:m3) = sense*(c(cixLow) - Prob.c_L(cixLow));
         g(m3+1:m) = sense*(Prob.c_U(cixUpp) - c(cixUpp));
      end
   end
end

% MODIFICATION LOG:
%
% 030427  hkh  Written
% 040117  med  fixed 2 bugs
% 060904  ango Avoid nlp_c call if only linearly constrained