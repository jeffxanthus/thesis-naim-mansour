% nlp_dcF is a utility routine when no special treatment of linear constraints
% are possible.
%
% nlp_dcF puts the derivative of all linear and nonlinear constraints 
% in one array.
%
% The order is:
% Linear    equalities
% Nonlinear equalities
% Linear    inequalities
% Nonlinear inequalities
%
% function dc=nlp_dcF(x, Prob, varargin)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Mar 19, 1998.  Last modified Dec 4, 2001.

function dc=nlp_dcF(x, Prob, varargin)

dcNL=nlp_dc(x, Prob, varargin{:}); 

dc=[];
if ~isempty(Prob.b_L)
   ix=find((Prob.b_L==Prob.b_U) & (abs(Prob.b_L) ~= Inf));
   if ~isempty(ix)
      dc=Prob.A(ix,:);
   end
end
if ~isempty(Prob.c_L)
   ix=find((Prob.c_L==Prob.c_U) & (abs(Prob.c_L) ~= Inf));
   if ~isempty(ix)
      dc=[dc;dcNL(ix,:)];
   end
end
if ~isempty(Prob.b_L)
   ix=find((Prob.b_L~=Prob.b_U) & (Prob.b_L ~= -Inf));
   if ~isempty(ix)
      dc=[dc;Prob.A(ix,:)];
   end
   ix=find((Prob.b_L~=Prob.b_U) & (Prob.b_U ~= Inf));
   if ~isempty(ix)
      dc=[dc; -Prob.A(ix,:)];
   end
end
if ~isempty(Prob.c_L)
   ix=find((Prob.c_L~=Prob.c_U) & (Prob.c_L ~= -Inf));
   if ~isempty(ix)
      dc=[dc; dcNL(ix,:)];
   end
   ix=find((Prob.c_L~=Prob.c_U) & (Prob.c_U ~= Inf));
   if ~isempty(ix)
      dc=[dc; -dcNL(ix,:)];
   end
end

% MODIFICATION LOG:
%
% 990626  hkh  Avoid feval
% 990909  hkh  Avoid computing structural information
% 001218  hkh  Change to use rows for constraints, columns for variables
% 011204  hkh  Bug in the handling of rows for constraints