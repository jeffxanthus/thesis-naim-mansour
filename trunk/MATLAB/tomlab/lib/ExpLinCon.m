% function Prob = ExpLinCon(Prob, Tol, Upper)
%
% Small code to generate linear constraints for ExpFit problem
%
% Should be called at the end of a zzz_prob routine for the
% exponential fitting problem, e.g. exp_prob.m
%
% Could also be called separately if using the TQ format
%
% The number of exponentials terms, Prob.ExpFit.p should be set beforehand
% It works independent of the separable/nonseparable strategy (Prob.LS.SepAlg) 
%
% Input parameters:
%
% Tol    Minimal tolerance between different lambda parameters, default 1E-7
% Upper  Upper bound for the difference between two lambda parameters, 
%        default 30. Not critical, should be large enough.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Oct 28, 1998. Last modified Jan 25, 2004.

function Prob = ExpLinCon(Prob, Tol, Upper)

if nargin < 2
   Tol = [];
   if nargin < 3
      Upper = [];
   end
end

if isempty(Tol),   Tol   = 1E-7; end
if isempty(Upper), Upper = 30; end

p      = Prob.ExpFit.p;
n      = length(Prob.x_0);
Upper  = 30;

if p > 1
   %if 0 % As possible linear constraints
   %   m=p*(p-1)/2;
   %   A=zeros(m,n);
   %   k=0;
   %   for i=1:p-1
   %       for j=i+1:p
   %           k=k+1;
   %           A(k,i)=-1;
   %           A(k,j)=1;
   %       end
   %   end
   %else
      % Simpler way, only p-1 inequalities
      m=p-1;
      A=zeros(m,n);
      k=0;
      for i=1:p-1
          A(i,i)=-1;
          A(i,i+1)=1;
      end
   %end
   if ~isempty(Prob.A)
      ix = find(Prob.b_L ~= Tol | Prob.b_U ~= Upper);
   else
      ix = [];
   end
   if ~isempty(ix)
      Prob.A   = [Prob.A(ix,:);A];
      Prob.b_L = [Prob.b_L(ix);Tol*ones(m,1)];
      Prob.b_U = [Prob.b_U(ix);Upper*ones(m,1)];
   else
      Prob.A   = A;
      Prob.b_L = Tol*ones(m,1);
      Prob.b_U = Upper*ones(m,1);
   end
else
   if ~isempty(Prob.A)
      ix = find(Prob.b_L ~= Tol | Prob.b_U ~= Upper);
      if ~isempty(ix)
         Prob.A   = Prob.A(ix,:);
         Prob.b_L = Prob.b_L(ix);
         Prob.b_U = Prob.b_U(ix);
      else
         Prob.A   = [];
         Prob.b_L = [];
         Prob.b_U = [];
      end
   end
end
Prob.mLin = size(Prob.A,1);

% MODIFICATION LOG
%
% 020422  hkh  Do not set Prob.N! Set in expInit
% 040125  hkh  Define field mLin