% nlp_cF is a utility routine when no special treatment of linear constraints
% are possible.
% nlp_cF puts all linear and nonlinear constraints in one array c.
%
% The order is:
% Linear    equalities
% Nonlinear equalities
% Linear    inequalities
% Nonlinear inequalities
%
% function [c, Pdim] = nlp_cF(x, Prob, varargin)
%
% In Pdim the different constraint dimensions are stored
%
% (1) Total number of constraints (not bounds) (=mcon)
% (2) Number of linear and nonlinear equalities (=me)
% (3) Number of linear equalities (=meq)
% (4) Number of linear inequalities (=mleq)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Mar 19, 1998.  Last modified Sep 9, 1999.

function [c, Pdim] = nlp_cF(x, Prob, varargin)

cNL=nlp_c(x, Prob, varargin{:});
cNL=cNL(:);

% Compute linear constraints
if ~isempty(Prob.A)
   Ax=Prob.A*x;
else
   Ax=[];
end

c=[];
if ~isempty(Prob.b_L)
   if isempty(Prob.b_U)
      ix=[];
   else
      ix=find((Prob.b_L==Prob.b_U) & (abs(Prob.b_L) ~= Inf));
   end
   if ~isempty(ix)
      c=[c;Ax(ix)-Prob.b_L(ix)];
   end
end

Pdim(3)=length(c);

if ~isempty(Prob.c_L)
   if isempty(Prob.c_U)
      ix=[];
   else
      ix=find((Prob.c_L==Prob.c_U) & (abs(Prob.c_L) ~= Inf));
   end
   if ~isempty(ix)
      c=[c;cNL(ix)-Prob.c_L(ix);];
   end
end

Pdim(2)=length(c);

if ~isempty(Prob.b_L)
   if isempty(Prob.b_U)
      ix=find(~isinf(Prob.b_L));
   else
      ix=find((Prob.b_L~=Prob.b_U) & (Prob.b_L ~= -Inf));
   end
   if ~isempty(ix)
      c=[c;Ax(ix)-Prob.b_L(ix)];
   end
   if ~isempty(Prob.b_U)
      ix=find((Prob.b_L~=Prob.b_U) & (Prob.b_U ~= Inf));
      if ~isempty(ix)
         c=[c; Prob.b_U(ix)-Ax(ix)];
      end
   end
elseif ~isempty(Prob.b_U)
   ix=find(~isinf(Prob.b_U));
   if ~isempty(ix)
      c=[c; Prob.b_U(ix)-Ax(ix)];
   end
end

Pdim(4)=length(c)-Pdim(2);

if ~isempty(Prob.c_L)
   if isempty(Prob.c_U)
      ix=find(~isinf(Prob.c_L));
   else
      ix=find((Prob.c_L~=Prob.c_U) & (Prob.c_L ~= -Inf));
   end
   if ~isempty(ix)
      c=[c; cNL(ix)-Prob.c_L(ix)];
   end
   if ~isempty(Prob.c_U)
      ix=find((Prob.c_L~=Prob.c_U) & (Prob.c_U ~= Inf));
      if ~isempty(ix)
         c=[c; Prob.c_U(ix)-cNL(ix)];
      end
   end
elseif ~isempty(Prob.c_U)
   ix=find(~isinf(Prob.c_U));
   if ~isempty(ix)
      c=[c; Prob.c_U(ix)-cNL(ix)];
   end
end

Pdim(1)=length(c);

% MODIFICATION LOG:
%
% 980826 hkh  Error in comments, Pdim 2nd output parameter.
%             Made cNL column vector.
% 981129 hkh  Safeguard if b_U or c_U empty
% 990626 hkh  Avoid feval