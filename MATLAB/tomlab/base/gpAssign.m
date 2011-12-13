% gpAssign is a direct way of setting up a Geometric Programming (GP) problem
% in the TOMLAB (TQ) format.
%
% The information is put into the TOMLAB input problem structure Prob.
%
% Prob = gpAssign(...)
%
% It is then possible to solve the problem using the colp_gp solver
% included in TOMLAB /GP. See the manual for more detailed information.
%
% -----------------------------------------------------
%
% GP minimization problem (primal):
%
%
%        min    g0(t)
%         t
%        s/t   gk(t) <= 1, for k = 1,2,...,p
%              t(i) > 0,   for i = 1,2,...,m
%
% where
%
%        g0(t) = sum{j=1,n(0)} [c(j)*t(1)^a(j,1)...t(m)^a(j,m)]
%
%        gk(t) = sum{j=n(k-1)+1,n(k)} [c(j)*t(1)^a(j,1)...t(m)^a(j,m)]
%        for k = 1,2,...,p
%
%
% The copl_gp program solves posynomial GP problems. The dual solved is:
%
%        max    prod( (c(j)/x(j))^x(j)) * prod(lambda(k)^lambda(k))
%         x
%        s/t    sum(x(i)) = 1,       for i=1,2,...,n(0)
%               sum(x(j)*a(j,i)) = 0 for i=1,2,...,m
%               x(j) >= 0, j=1,2,...n(p)
%
% where lambda(k) = sum(x(j),j=n(k-1)+1:n(k)), k=1,2,...p
%
%
% -----------------------------------------------------
%
% Syntax of gpAssign:
%
% function Prob = gpAssign(nterm, coef, A, Name)
%
% INPUT (One parameter c must always be given)
%
% nterm:       Number of product terms in the objective followed by the
%              number of product terms in the constraints.
%              Example: nterm = [2;3;3;6] means 2 product terms in the
%              objective and 3 in the first and second constraint. The
%              third constraint has 6 terms. In the problem
%              description, nterm is called n.
%
% coef:        Positive coefficients c(j), j=1,...,n(p).
% A:           Matrix of exponents in primal problem. Transposed
%              linear constraints matrix in dual problem.
%
% Name         Problem name.

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written May 9, 2005. Last modified Oct 23, 2009.

function Prob = gpAssign(nterm, c, A, Name)

if nargin < 4
   Name=[];
   if nargin < 3
      A=[];
      if nargin < 2
         c=[];
         if nargin < 1
            nterm=0;
end, end, end, end

global MAX_x % Max number of variables/constraints to print

% CHECKS

if size(A,1) ~= length(c)
   fprintf('A must have the same number of rows\n');
   error('as c has entries.');
end

if ~all(c>=0)
   fprintf('c can only have positive entries.\n');
   error('Chech the inputs and try to convert the problem.');
end

A = sparse(A);

probType = checkType('gp');

Prob=ProbDef(1);

Prob.GP.nterm = nterm(:);
Prob.GP.coef = c(:);
Prob.GP.A    = A;

Prob.GP.b_U    = [zeros(size(A,2),1)];
Prob.GP.b_L    = Prob.GP.b_U;

Prob.P      = 1;
n           = size(A,2);
Prob.N      = n;
Prob.x_L    = zeros(n,1);
Prob.x_U    = inf*ones(n,1);

Prob.probType = probType;
Prob.probFile = 0;

Prob.mLin     = 0;

Prob.c_L     = -inf*ones(length(Prob.GP.nterm)-1,1);
Prob.c_U     = ones(length(Prob.GP.nterm)-1,1);
Prob.mNonLin = length(Prob.c_U);

% Set Print Level to 0 as default
Prob.PriLevOpt=0;

if isempty(MAX_x)
   MAX_x=20;
end

Prob = tomFiles(Prob, 'gp_f', [], [], 'gp_c');
Prob.Name = deblank(Name);

% MODIFICATION LOG
%
% 050509  med  Created
% 050519  med  Modified to new format
% 050520  med  Modified to handle transposed A
% 050526  med  Updated size of problem, Prob.N
% 050527  frhe Updated problem description and help text
% 050608  med  Updated so general solver can run
% 070627  ango Help text updated
% 091023  med  Help updated