% function qpProb = lls2qp(Prob, IntVars)
%
% Converts an lls problem to a new problem based on the formula below.
% Only the objective function is affected.
% The problem can be of any type with an LLS objective.
%
%     0.5 * |y - Cx| = 0.5 * (y - Cx)'*(y - Cx) =
%     = 0.5 * (y'y - 2y'Cx + x'C'Cx) =>
%     => 0.5 * y'y + 0.5 * x'Fx + c'x for
%
%     F = C'C and c' = -y'C, the problem also has a constant part 0.5 * y'y
%
% REQUIRED INPUTS:
%
% Prob.LS.C      The linear matrix in 0.5 * |y - Cx|
% Prob.LS.y      The constant vector in 0.5 * |y - Cx|
%
% OPTIONAL:
%   IntVars:
%             If empty, Check Prob.MIP.IntVars
%             If islogical(IntVars) (=all elements are 0/1), then
%             1 = integer variable, 0 = continuous variable.
%             If any element >1, IntVars is the indices for integer variables
%
% OUTPUTS:
%
% qpProb         The converted problem
%
% If the problem is a linear least squares problem a qp problem is created.
% The new problem may have integer variables. Create the problem with
% llsAssign then use this routine.
%
% If the problem has nonlinear constraints an nlp is created.
% The new problem may have integer variables. Create the problem with
% conAssign or minlpAssign, then set the fields Prob.LS.C and Prob.LS.y

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Feb 21, 2005.   Last modified Jul 17, 2009.

function qpProb = lls2qp(Prob, IntVars)

if nargin < 2
   IntVars = [];
   if nargin < 1
      error('The input Prob is required');
   end
end

% Check that everything exist
if isempty(Prob.LS)
   error('Field Prob.LS is empty');
else
   if isfield(Prob.LS, 'C')
      if isempty(Prob.LS.C)
         error('Field Prob.LS.C is empty');
      end
   end
   if isfield(Prob.LS, 'y')
      if isempty(Prob.LS.y)
         error('Field Prob.LS.y is empty');
      end
   end
end

% Copy problem
qpProb = Prob;

if isempty(IntVars)
   % IntVars may be in Prob.MIP.IntVars already
   % Integer variables
   IntVars  = DefPar(Prob.MIP,'IntVars',[]);
end

n  = Prob.N;
% Logical vector for integers
IV = false(n,1);

if isempty(IntVars)
   % No binary variables B or integer variables of type I
elseif any(IntVars==0) | all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > n)
      error('arbf: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end
IntVars = find(IV);
qpProb.MIP.IntVars = IntVars(:);

% Set correct problem type
if Prob.probType == checkType('lls')
   if isempty(IntVars)
      qpProb.probType = checkType('qp');
   else
      qpProb.probType = checkType('miqp');
   end
else
   if isempty(IntVars)
      qpProb.probType = checkType('con');
   else
      qpProb.probType = checkType('minlp');
   end
end

% Check sizes on inputs
if size(Prob.LS.C, 2) ~= Prob.N
   error('The number of columns in Prob.LS.C does not match problem size');
end

Prob.LS.y = Prob.LS.y(:);

if length(Prob.LS.y) ~= size(Prob.LS.C, 1)
   error('The length of Prob.LS.y does not match problem size');
end

% Remove all traces of an LLS problem
qpProb.LS = [];

% Add identifier to the name telling the problem has been converted
qpProb.Name = [qpProb.Name ' (lls2qp)'];

qpProb.QP.F = Prob.LS.C'*Prob.LS.C;
qpProb.QP.c = [-Prob.LS.y'*Prob.LS.C]';
qpProb.fConstant = 0.5*(Prob.LS.y'*Prob.LS.y);

qpProb.FUNCS.f = 'qp_f';
qpProb.FUNCS.g = 'qp_g';
qpProb.FUNCS.H = 'qp_H';

% MODIFICATION LOG:
%
% 050221 frhe File written
% 050422 frhe fConstant added
% 051216 med  Function released and polished
% 070222 hkh  Revise IntVars handling (code not OK), use new format
% 070223 hkh  Set n=Prob.N for use in IntVars and IV definition
% 090717 med  fConstant calculation updated
