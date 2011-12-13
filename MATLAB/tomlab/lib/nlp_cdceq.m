% nlp_cdceq.m
%
% function [c, dc]=nlp_cdceq(x, Prob, varargin)
%
% nlp_cdceq calls the TOMLAB gateway function nlp_c, 
% which evaluates the constraints at x for the test problem P (Prob.P).
%
% if nargout > 2, it calls the TOMLAB gateway function nlp_dc, 
% which computes the gradient for all constraints at x for the test problem P
%
% The output is divided into inequality and equality constraints
% and transformed to fit the formulation:
%
%  c(x) <= 0, c(x)== 0
%
% nlp_cdceq is used when calling OPT TB 2.0

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Mar 5, 1997.   Last modified April 27, 2003.

function [c, ceq, dc, dceq] = nlp_cdceq(x, Prob, varargin)

c= nlp_c( x, Prob, varargin{:});

if ~isempty(Prob.c_L)
   ixE=find(Prob.c_L==Prob.c_U & ~isinf(Prob.c_L) );
   if ~isempty(ixE)
      ceq=c(ixE);
   else
      ceq=[];
   end
else
   ceq=[];
   ixE=[];
end
if ~isempty(Prob.c_L)
   ixL=find(Prob.c_L~=Prob.c_U & ~isinf(Prob.c_L));
   ixU=find(Prob.c_L~=Prob.c_U & ~isinf(Prob.c_U));
   if isempty(ixL)
      c=c(ixU)-Prob.c_U(ixU);
   else
      if isempty(ixU)
         c=Prob.c_L(ixL)-c(ixL);
      else
         c=[Prob.c_L(ixL)-c(ixL);c(ixU)-Prob.c_U(ixU)];
      end
   end
else
   ixL=[]; ixU=[];
end

if nargout > 2
   dc=(nlp_dc(x, Prob, varargin{:}));  
   if ~isempty(ixE)
      dceq=dc(ixE,:)';
   else
      dceq=[];
   end
   if isempty(ixL)
      dc=dc(ixU,:)';
   else
      if isempty(ixU)
         dc=-dc(ixL,:)';
      else
         dc=[-dc(ixL,:); dc(ixU,:)]';
      end
   end
end

% opt tbx syntax: Each column in dc now corresponds to a constraint.

% MODIFICATION LOG:
%
% 990626  hkh  Avoid feval
% 990705  hkh  Make dc a full matrix
% 030427  hkh  Change tests on inf to use isinf