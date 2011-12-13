% mipFeasible.m
%
% mipFeasible finds feasible initial grid, w.r.t. linear constraints
% 
% It tries to solve LP/MILP problems with linear different objectives
%
% function X = mipFeasible(Prob, M, nSample, Method, PriLev)
%
% INPUT:
%
% Prob      Problem structure
%           If Prob.CGO.X is defined, generate nSample additional points
%
% M         Number of trial points to generate
%
% nSample   Try to find nSample feasible points
%
% Method    Method. Now only Method = 1 implemented
%
% PriLev    PriLev = 0, no output, =1 summary, =2 Detailed output
%
% OUTPUT:
% X         Sample point matrix, d x M

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written June 14, 2008.   Last modified Sep 2, 2009.

function X = mipFeasible(Prob, M, nSample, Method, PriLev)

if nargin < 5 
   PriLev = []; 
   if nargin < 4 
      Method = []; 
   end
end

if isempty(Method), Method = 1; end
if isempty(PriLev), PriLev = 0; end

x_L     = Prob.x_L;
x_U     = Prob.x_U;
d       = Prob.N;
IntVars = Prob.MIP.IntVars;

% Generate M objective functions
if isfield(Prob.CGO,'X')
   X       = Prob.CGO.X;
else
   X       = [];
end

% Genereate nFac as many, to assure full rank
if any(Prob.b_L==Prob.b_U)
   nFac    = 4;
else
   nFac    = 2;
end
nSample = nSample*nFac;

if isempty(X)
   X       = zeros(d,nSample);
   k       = 0;
else
   k       = size(X,2);
   X       = [X,zeros(d,nSample)];
   nSample = nSample+k;
end
% Avoid use of cplex
%if isempty(IntVars)
%   SolverMIP = GetSolver('lp');
%else
%   SolverMIP = GetSolver('mip');
%end
% Should use a solver in Base module
SolverMIP='milpSolve';

%SolverMIP='xa';
%SolverMIP='bpqd';
if PriLev > 0
   fprintf('mipFeasible: Find feasible initial points using %s\n',SolverMIP);
   fprintf('           : Method = %d. Trial points M = %d\n',Method,M);
end

if Method == 1
   if Prob.mLin > 0
      ixEQ = find(Prob.b_L == Prob.b_U);
      % ixEQ = [1:Prob.mLin];
      nEQ = length(ixEQ);
      if nEQ > 0
         [Q R] = qr(Prob.A(ixEQ,:)');
	 pRank = rank(R);
         nNS   = size(Q,2)-pRank;
         % Generate nNS x vectors
         xMP   =(x_U-x_L)/2;
         V = xMP*ones(1,nNS) + xMP*ones(1,nNS).*Q(:,pRank+1:pRank+nNS);
         % Z = V;
         V(IntVars,:) = round(V(IntVars,:));
         %for i=1:nNS
         %    D(i) = norm(V(:,i)-Z(:,i),2);
         %end
         %xprint(D,'D:')
         X(:,k+1:k+nNS) = V;
         k = k+nNS;
         if PriLev > 1
            fprintf('Use QR-NullSpace, Found %2d - nSample %d\n',nNS,nSample);
         end
         % NOT X = [X, V(:,nEQ+1:end)];
      end
   else
      % nEQ = 0; % NOT USED
      
      [X,nCon] = randomDesigns(3,x_L,x_U,1,round(nSample-k)/4,[],X(:,1:k),...
                                             M,Prob,0,0,100+PriLev);
      return
   end
   VW = DefPar(Prob.MIP,'VarWeight',[]);
   ProbMIP = mipAssign(zeros(d,1), Prob.A, Prob.b_L, Prob.b_U, x_L, x_U, ...
                Prob.x_L,'mipFeasible',[],[],IntVars, VW);
   VW = VW(:);
   ProbMIP.PriLevOpt = 0;
   ProbMIP = ProbCheck(ProbMIP,SolverMIP,ProbMIP.probType);
   ProbMIP.optParam.IterPrint = 0;

   if strcmpi('cplex',SolverMIP)

      if ~isempty(VW)
         [i1,i2]=sort(VW);
         Priority(i2(i2)) = -i2;
         Priority = Priority - min(Priority);
         ProbMIP.CPLEX.BranchPrio = Priority;
         xprinti(Priority,'Prio:')
      end
   end
   if 1 & ~isempty(VW)
      ProbMIP.QP.c    = VW;
      R               = tomRunFast(SolverMIP,ProbMIP);
      [X,k]           = addNewPoint(R,X,k);
      if PriLev > 1
         fprintf('k=%3d. Prio %12.6f.\n',k,sum(X(:,k).*VW));
         fprintf('UseVW   0,        Found %2d - nSample %d\n',k,nSample);
      end
   end
   if 0 & ~isempty(VW)
      mx = 10*max(abs(VW));
      for i=1:d 
          ProbMIP.QP.c(i) = VW(i)+mx; 
          R               = tomRunFast(SolverMIP,ProbMIP);
          [X,k]           = addNewPoint(R,X,k);
          fprintf('k=%3d. Prio %12.6f.\n',k,sum(X(:,k).*VW));
          ProbMIP.QP.c(i) = VW(i); 
          if k >= nSample, break; end
      end
      if PriLev > 1
         fprintf('UseVW   1,        Found %2d - nSample %d\n',k,nSample);
      end
      if k >= nSample, return; end
      % m1 = mean(VW); % NOT USED
      s1 = std(VW);
      for i=1:0 
          ProbMIP.QP.c    = VW+s1*randn(Prob.N,1); 
          R               = tomRunFast(SolverMIP,ProbMIP);
          [X,k]           = addNewPoint(R,X,k);
          ProbMIP.QP.c    = VW;
          if k >= nSample, break; end
      end
      if PriLev > 1
         fprintf('UseVW std,        Found %2d - nSample %d\n',k,nSample);
      end
      if k >= nSample, return; end
   end

   for i=1:d 
       ProbMIP.QP.c(i) = 1; 
       R               = tomRunFast(SolverMIP,ProbMIP);
       k0              = k;
       [X,k]           = addNewPoint(R,X,k);
       if PriLev > 1 & k > k0 & ~isempty(VW)
          fprintf('k=%3d. Prio %12.6f.\n',k,sum(X(:,k).*VW));
       end
       ProbMIP.QP.c(i) = 0; 
       if k >= nSample, break; end
   end
   if PriLev > 1
      fprintf('Usevar  1,        Found %2d - nSample %d\n',k,nSample);
   end
   if k >= nSample, return; end
   for i=1:d 
       ProbMIP.QP.c(i) = -1; 
       R               = tomRunFast(SolverMIP,ProbMIP);
       k0              = k;
       [X,k]           = addNewPoint(R,X,k);
       if PriLev > 1 & k > k0 & ~isempty(VW)
          fprintf('k=%3d. Prio %12.6f.\n',k,sum(X(:,k).*VW));
       end
       ProbMIP.QP.c(i) = 0; 
       if k >= nSample, break; end
   end
   if PriLev > 1
      fprintf('Usevar -1,        Found %2d - nSample %d\n',k,nSample);
   end
   if k >= nSample, return; end

   for i=1:d 
       ProbMIP.QP.c(i) = 1; 
       for j=i+1:d 
           ProbMIP.QP.c(j) = 1; 
           R               = tomRunFast(SolverMIP,ProbMIP);
           k0              = k;
           [X,k]           = addNewPoint(R,X,k);
           if PriLev > 1 & k > k0 & ~isempty(VW)
              fprintf('k=%3d. Prio %12.6f.\n',k,sum(X(:,k).*VW));
           end
           ProbMIP.QP.c(j) = 0; 
           if k >= nSample, break; end
       end
       ProbMIP.QP.c(i) = 0; 
       if k >= nSample, break; end
       for j=i+1:d 
           ProbMIP.QP.c(j) = -1; 
           R               = tomRunFast(SolverMIP,ProbMIP);
           k0              = k;
           [X,k]           = addNewPoint(R,X,k);
           if PriLev > 1 & k > k0 & ~isempty(VW)
              fprintf('k=%3d. Prio %12.6f.\n',k,sum(X(:,k).*VW));
           end
           ProbMIP.QP.c(j) = 0; 
           if k >= nSample, break; end
       end
       ProbMIP.QP.c(i) = 0; 
       if k >= nSample, break; end
   end
   if PriLev > 1
      fprintf('Usevar  2,        Found %2d - nSample %d\n',k,nSample);
   end
   if k >= nSample, return; end

   for i=1:d 
       ProbMIP.QP.c(i) = -1; 
       for j=i+1:d 
           ProbMIP.QP.c(j) = 1; 
           R               = tomRunFast(SolverMIP,ProbMIP);
           k0              = k;
           [X,k]           = addNewPoint(R,X,k);
           if PriLev > 1 & k > k0 & ~isempty(VW)
              fprintf('k=%3d. Prio %12.6f.\n',k,sum(X(:,k).*VW));
           end
           ProbMIP.QP.c(j) = 0; 
           if k >= nSample, break; end
       end
       ProbMIP.QP.c(i) = 0; 
       if k >= nSample, break; end
       for j=i+1:d 
           ProbMIP.QP.c(j) = -1; 
           R               = tomRunFast(SolverMIP,ProbMIP);
           k0              = k;
           [X,k]           = addNewPoint(R,X,k);
           if PriLev > 1 & k > k0 & ~isempty(VW)
              fprintf('k=%3d. Prio %12.6f.\n',k,sum(X(:,k).*VW));
           end
           ProbMIP.QP.c(j) = 0; 
           if k >= nSample, break; end
       end
       ProbMIP.QP.c(i) = 0; 
       if k >= nSample, break; end
   end
   if PriLev > 1
      fprintf('Usevar -2,        Found %2d - nSample %d\n',k,nSample);
   end
   if k >= nSample, return; end

   ProbMIP.QP.c = ones(d,1); 
   R            = tomRunFast(SolverMIP,ProbMIP);
   k0           = k;
   [X,k]        = addNewPoint(R,X,k);
   if PriLev > 1 & k > k0 & ~isempty(VW)
      fprintf('k=%3d. Prio %12.6f.\n',k,sum(X(:,k).*VW));
   end
   if PriLev > 1
      fprintf('Usevar %2d,        Found %2d - nSample %d\n',d,k,nSample);
   end
   if k >= nSample, return; end

   ProbMIP.QP.c = -ones(d,1); 
   R            = tomRunFast(SolverMIP,ProbMIP);
   k0           = k;
   [X,k]        = addNewPoint(R,X,k);
   if PriLev > 1 & k > k0 & ~isempty(VW)
      fprintf('k=%3d. Prio %12.6f.\n',k,sum(X(:,k).*VW));
   end
   if PriLev > 1
      fprintf('Usevar %2d,        Found %2d - nSample %d\n',-d,k,nSample);
   end
   if k >= nSample, return; end

   ProbMIP.QP.c = zeros(d,1); 

   if d > 2
      p = ceil(M/(d-2));
      for j=1:p
      for i=d:-1:3
          [v,ix] = sort(rand(i,1));
          ProbMIP.QP.c(ix(1:i)) = 10*rand(i,1); 
          R                     = tomRunFast(SolverMIP,ProbMIP);
          k0                    = k;
          [X,k]                 = addNewPoint(R,X,k);
          if PriLev > 1 & k > k0 & ~isempty(VW)
             fprintf('k=%3d. Prio %12.6f.\n',k,sum(X(:,k).*VW));
          end
          ProbMIP.QP.c(ix(1:i)) = 0;
          if k >= nSample, break; end
      end
      if PriLev > 1
         fprintf('Cycle  %2d of %3d, Found %2d - nSample %d\n',j,p,k,nSample);
      end
      if k >= nSample, break; end
      end
      if k >= nSample, return; end
   end

   for i=1:M
       ProbMIP.QP.c = rand(d,1)-0.5; 
       R            = tomRunFast(SolverMIP,ProbMIP);
       k0           = k;
       [X,k]        = addNewPoint(R,X,k);
       if PriLev > 1 & k > k0 & ~isempty(VW)
          fprintf('k=%3d. Prio %12.6f.\n',k,sum(X(:,k).*VW));
       end
       if k >= nSample, break; end
   end
   if PriLev > 0
      fprintf('           : Found k = %2d points - needed %d\n',k,nSample/nFac);
   end
end

% if k < nSample , return only the unique points found
X = X(:,1:k);


function [X,k] = addNewPoint(R,X,k)

if R.ExitFlag ~= 4
   x_k = R.x_k;
   for i = 1:size(x_k,2)
       if k== 0 | tomsol(30,x_k(:,i),X(:,1:k)) > 1E-7
          k      = k+1;
          X(:,k) = x_k(:,1);
       %else
          %fprintf('k %d ',k);
          %xprinti(x_k(:,1),'x:');
       end
   end
end

% MODIFICATION LOG:
%
% 080614  hkh  Written
% 080615  hkh  Add PriLev
% 080711  hkh  Use Prob.MIP.VarWeight, if defined. 
% 080711  hkh  Redefine VarWeight for CPLEX, defined in Prob.CPLEX.BranchPrio
% 080716  hkh  Add negative coefficients as well. Only return the k feasible X
% 090429  hkh  if no eq, call randomDesigns to get new trial points
% 090824  hkh  Minor mlint revision, avoid some statements
% 090901  hkh  PriLev > 1 more most output, only summary with PriLev > 0
% 090902  hkh  nFac=2 (nFac*nSample) if no linear equalities, use milpSolve
% 090902  hkh  Correction of QR Null space part
