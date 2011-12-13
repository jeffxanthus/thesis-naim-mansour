% nlp2_cdceq.m
%
% function [c, ceq, dc, dceq] = nlp2_cdceq(x, Prob, varargin)
%
% nlp2_cdceq returns the constraints c(x), and the constraint gradients dc(x)
% for a nonlinear problem.
%
% nlp2_cdceq calls the user function that returns the constraint, and possibly
% the user function that returns the constraint gradient vectors
%
% if nargout > 1, nlp_cdceq also computes the constraint gradients
%
% nlp2_cdceq is used when implementing the OPTIM TB 2.x compatibility interface
%
% The output is divided into inequality and equality constraints
% and transformed to fit the formulation:
%
%  c(x) <= 0, c(x)== 0

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written July 6, 1999.   Last modified Aug 14, 2006.

function [c, ceq, dc, dceq] = nlp2_cdceq(x, Prob, varargin)

%nargin;
global n_c n_dc NARG

cFunc=Prob.FUNCSX.c;

if isempty(cFunc)
   c=[]; ceq=[]; dc=[]; dceq=[];
   return
end

n_c = n_c + 1;

if NARG(4) > 1
   c = feval(cFunc, x, Prob, varargin{:});
%elseif xnargin(cFunc) > 1
%   c = feval(cFunc, x, Prob);
else
   c = feval(cFunc, x);
end

% Split c into two parts

if ~isempty(Prob.c_L)
   ixE=find((Prob.c_L==Prob.c_U) & (abs(Prob.c_L) ~= Inf));
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
   ixL=find((Prob.c_L~=Prob.c_U) & (Prob.c_L ~= -Inf));
   ixU=find((Prob.c_L~=Prob.c_U) & (Prob.c_U ~= Inf));
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
   dcFunc=Prob.FUNCSX.dc;
   if ~isempty(dcFunc)
      n_dc = n_dc + 1;
      if NARG(5) > 1
         dc = (feval(dcFunc, x, Prob, varargin{:}));
      %elseif xnargin(dcFunc) > 1
      %   dc = (feval(dcFunc, x, Prob));
      else
         dc = (feval(dcFunc, x));
      end
      % REMARK. Transposed matrix!!! 
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
   else
      dc   = [];
      dceq = [];
   end

end

% opt tbx syntax: Each column in dc now corresponds to a constraint.

% MODIFICATION LOG:
%
% 990626  hkh  Avoid feval
% 990705  hkh  Make dc a full matrix
% 000831  hkh  Avoid infinite recursion, change logic
% 030129  hkh  Major revision for v4.0
% 040125  hkh  Transpose missing for dceq=dc(ixE,:)';, line 81
% 060814  med  FUNCSX used for callbacks instead