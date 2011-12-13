% LagMult computes Lagrange multiplier estimates.
%
% If requested, also computes the reduced (projected) gradient
%
% function [v_k, Zv, P, Z_L, cErr, ceq, cineq, gProj] = LagMult(Prob, Result);
%
% INPUT:
% Prob        The TOMLAB problem input structure
%             Using fields Prob.x_L; Prob.x_U; Prob.A; Prob.b_L; Prob.b_U;
%             Prob.c_L; Prob.c_U; Prob.optParam.eps_Rank;
%             Prob.optParam.xTol; Prob.optParam.bTol; Prob.optParam.cTol;
% Result      The TOMLAB optimization result output structure
%             Using fields: Result.x_k; Result.g_k; Result.c_k;
%
% OUTPUT:
% v_k         The full Lagrange multiplier vector
%             Order: simple bounds, linear constraints, nonlinear constraints
% Zv          Description of all active constraints
% P           Status indicator for each constraint
%             0 = inactive, -1 = active on lower bound
%             2 = equality,  1 = active on upper bound
% Z_L         The matrix with derivatives of the active constraints
% cErr        Relative error in nonlinear constraints
% ceq         Number of too large errors in nonlinear equality constraints
% cineq       Number of too large errors in nonlinear inequality constraints
% gProj       Reduced (Projected) gradient (I-Q*Q')g_k, where Z_L=Q*R*E'

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.1.0$
% Written Dec 8, 1998.      Last modified April 17, 2008.

function [v_k, Zv, P, Z_L, cErr, ceq, cineq, gProj] = LagMult(Prob, Result)

if nargin < 2
   error('LagMult needs two inputs, structures Prob and Result');
end

Zv    = [];

x_k  = Result.x_k;
g_k  = Result.g_k;
c_k  = Result.c_k;
[n,l]= size(x_k);
if n ~= Prob.N
   fprintf('\n\n');
   fprintf('Number of variables, Prob.N                %d\n',Prob.N);
   fprintf('Number of computed nonlinear variables x_k %d\n',n);
   fprintf('\n');
   fprintf('WARNING!!! Possible input errors !!! \n');
   fprintf('\n');
end
if l > 1
   x_k = x_k(:,1);
end

Z_L  = zeros(n,0);
xTol = Prob.optParam.xTol;
bTol = Prob.optParam.bTol;
cTol = Prob.optParam.cTol;

if ~isempty(Prob.A)
   Ax = Prob.A*x_k;     
else
   Ax = zeros(0,1);
end

mA = Prob.mLin;
% mA=size(Prob.A,1);
m   = Prob.mNonLin;
[mNL,l] = size(c_k);
if l > 1
   c_k = c_k(:,1);
end
if m ~= mNL
   fprintf('\n\n');
   fprintf('Number of nonlinear constraints, Prob.mNonLin %d\n',m);
   fprintf('Number of computed nonlinear constraints      %d\n',mNL);
   fprintf('\n');
   fprintf('WARNING!!! Possible input errors !!! \n');
   fprintf('\n');
   m = mNL;
end

v_k = zeros(n+mA+m,1);
P   = zeros(n+mA+m,1);

x_L = [Prob.x_L;-inf*ones(n-length(Prob.x_L),1)];
x_U = [Prob.x_U; inf*ones(n-length(Prob.x_U),1)];
% Check simple bounds
for i=1:n
   if x_L(i)==x_U(i)
      z=zeros(n,1);
      z(i)=1;
      Z_L=[Z_L z];
      ss=sprintf('xVarEq  %d',i);
      Zv=str2mat(Zv,ss);
      v_k(i) = 1;
      P(i)=2;
   elseif norm(x_k(i)-x_L(i)) < xTol
      z=zeros(n,1);
      z(i)=1;
      Z_L=[Z_L z];
      ss=sprintf('xVarLow %d',i);
      Zv=str2mat(Zv,ss);
      v_k(i) = 1;
      P(i)=-1;
   elseif norm(x_k(i)-x_U(i)) < xTol
      z=zeros(n,1);
      z(i)=-1;
      Z_L=[Z_L z];
      ss=sprintf('xVarUpp %d',i);
      Zv=str2mat(Zv,ss);
      v_k(i) = 1;
      P(i)=1;
   end
end

% Check linear constraints
b_L = [Prob.b_L;-inf*ones(mA-length(Prob.b_L),1)];
b_U = [Prob.b_U; inf*ones(mA-length(Prob.b_U),1)];
if mA > 0
   for i=1:length(b_L)
      if b_L(i)==b_U(i)
         Z_L=[Z_L full(Prob.A(i,:))'];
         ss=sprintf('LinEq   %d',i);
         Zv=str2mat(Zv,ss);
         v_k(n+i) = 1;
         P(n+i)=2;
      elseif norm(Ax(i)-b_L(i)) < bTol
         Z_L=[Z_L full(Prob.A(i,:))'];
         ss=sprintf('LinLow  %d',i);
         Zv=str2mat(Zv,ss);
         v_k(n+i) = 1;
         P(n+i)=-1;
      elseif norm(Ax(i)-b_U(i)) < bTol
         Z_L=[Z_L -full(Prob.A(i,:))'];
         ss=sprintf('LinUpp  %d',i);
         Zv=str2mat(Zv,ss);
         v_k(n+i) = 1;
         P(n+i)=1;
      end
   end
end

% Check nonlinear constraints
if m > 0
   c_L = [Prob.c_L;-inf*ones(m-length(Prob.c_L),1)];
   c_U = [Prob.c_U; inf*ones(m-length(Prob.c_U),1)];
   dc_k = Result.cJac;
   if isempty(dc_k)
      JAC = 0;
   else
      JAC = 1;
   end
   cErr=zeros(m,1);
else
   cErr=[];
   JAC = 1;
end

ceq=0;
cineq=0;
for i=1:m
   if c_L(i)==c_U(i)
      cErr(i)=-abs(c_k(i)-c_L(i))/max(1,abs(c_L(i)));
      ceq=ceq+(abs(cErr(i)) > cTol);
      if JAC
         Z_L=[Z_L full(dc_k(i,:)')];
      end
      ss=sprintf('NonLEQ  %d',i);
      Zv=str2mat(Zv,ss);
      v_k(n+mA+i) = 1;
      P(n+mA+i)=2;
   elseif (c_k(i)-c_L(i)) < cTol*max(1,abs(c_L(i)));
      cErr(i)=(c_k(i)-c_L(i))/max(1,abs(c_L(i)));
      if JAC
         Z_L=[Z_L full(dc_k(i,:)')];
      end
      cineq=cineq+ (abs(cErr(i)) > cTol);
      ss=sprintf('NonLLow %d',i);
      Zv=str2mat(Zv,ss);
      v_k(n+mA+i) = 1;
      P(n+mA+i)=-1;
   elseif c_U(i)-c_k(i) < cTol*max(1,abs(c_U(i)));
      cErr(i)=(c_U(i)-c_k(i))/max(1,abs(c_L(i)));
      if JAC
         Z_L=[Z_L -full(dc_k(i,:)')];
      end
      cineq=cineq+ (abs(cErr(i)) > cTol);
      ss=sprintf('NonLUpp %d',i);
      Zv=str2mat(Zv,ss);
      v_k(n+mA+i) = 1;
      P(n+mA+i)=1;
   end
end
if ~JAC
   gProj=[];
   return
end

% Estimate Lagrange multipliers safely
% v = Z_L \ g_k

m2=min(size(Z_L));

%if 0 &  m2 < size(Z_L,2)
%   fprintf('TOMLAB found %d active constraints. ',size(Z_L,2))
%   fprintf('Estimate Lagrange multipliers for %d first\n',m2);
%end

if m2 > 0
   if m2 < size(Z_L,2)
      % In case of redundancy, put the linear and nonlinear constraints first
      nV=sum(P(1:n)~=0);
      jx=[nV+1:m2,nV:-1:1];
      jx=1:m2;
   else
      jx=1:m2;
   end
   [Q, R, E, pRank] = ComputeQR(Z_L(:,jx), Prob.optParam.eps_Rank);

   if isempty(g_k) 
      v_k=[];
      gProj=[];
   else
      if pRank==0
         v = zeros(length(jx),1);
      else
         v = tomsol(6, Q, R, E, pRank, g_k);
      end

      % Compute projected gradient, if requested

      if nargout > 7
         gProj = g_k - Q(:,1:pRank)*(Q(:,1:pRank)'*g_k);
      end

      ix=find(v_k);
      % Set the m2 multipliers into full v_k
      v_k(ix(jx)) = v;
   end
else
   gProj=g_k;
end

% MODIFICATION LOG
%
% 981210  hkh  Test if Prob.A is empty
% 990810  hkh  Switching from SVD to QR, handling large and sparse.
%              Computing projected gradient efficiently
% 990909  hkh  Do not call nlp_dc, use cJac in structure
% 001218  hkh  Change to transposed dc_k, each row one constraint
% 030309  hkh  Force efficient computation of projected gradient, add comments
% 031126  hkh  Check if Result.cJac [], avoid crash computing projected gradient
% 040102  hkh  Safe guard output for length differences in constraints
% 040102  hkh  Check for length differences between results and fields
% 040331  hkh  Check for multiple solutions x_k and c_k
% 080417  hkh  Change to relative constraint violation cErr instead of absolute
