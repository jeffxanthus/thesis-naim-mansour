% cpTransf.m
%
% Different transformations of nonlinear programs
%
% function [AA, bb, meq] = cpTransf(Prob,TransfType,makeEQ,LowInf)
%
% TransfType == 1. Transform convex program
%
%        min   f(x).  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U     Fixed    variables: Set x_L==x_U
%              b_L <= A x  <= b_U     Equality equations: Set b_L==b_U
%
%     to
%        min   f(x-x_L).  x in R^n
%         x
%        s/t    AA (x-x_L)  <= bb  (The first meq are equalities)
%                  (x-x_L)  >= 0
%
%     Translate back with (fixed variables do not change their values):
%        x(~x_L==x_U) = (x-x_L) + x_L(~x_L==x_U)
%
% TransfType == 2. Transform convex program
%
%        min   f(x).  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U     Fixed    variables: Set x_L==x_U
%              b_L <= A x  <= b_U     Equality equations: Set b_L==b_U
%
%     to
%        min   f(x).  x in R^n
%         x
%        s/t         AA x  <= bb  (The first meq are equalities)
%              x_L <=   x  <= x_U
%
% TransfType == 3. Transform convex program
%
%        min   f(x).  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U     Fixed    variables: Set x_L==x_U
%              b_L <= A x  <= b_U     Equality equations: Set b_L==b_U
%
%     to
%        min   f(x).  x in R^n
%         x
%        s/t         AA x  <= bb  (The first meq are equalities)
%              x_L <=   x
%
% INPUT PARAMETERS
% Fields in Prob:
%   QP.c:   The vector c in c'x
%   A:      The linear constraint matrix
%   b_L:    The lower bounds for the linear constraints
%   b_U:    The upper bounds for the linear constraints
%   x_L:    Lower bounds on x
%   x_U:    Upper bounds on x
%
%   TransfType  == 1,2 or 3
%   MakeEQ  Flag, if set make standard form (all equalities)
%   LowInf  Variables -Inf or variables < LowInf are set to LowInf before
%           transforming the problem.  Default = -1E4.
%           abs(LowInf) are limit if upper bound variables are to be used
%
% OUTPUT PARAMETERS
%   AA:     The expanded linear constraint matrix
%   bb:     The expanded upper bounds for the linear constraints
%   meq:    The first meq equations are equalities

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Nov 17, 1998.  Last modified Aug 13, 2009.

function [AA, bb, meq] = cpTransf(Prob,TransfType,makeEQ,LowInf)

if nargin < 4
   LowInf=[];
   if nargin < 3
      makeEQ=[];
      if nargin < 2
         TransfType=[];
      end
   end
end

%INPUT scalar: LowInf makeEQ TransfType

if isempty(LowInf),     LowInf=-1E4;  end
if isempty(makeEQ),     makeEQ=0;     end
if isempty(TransfType), TransfType=1; end

A=Prob.A;
[m n] = size(A);
b_L=Prob.b_L(:);
b_U=Prob.b_U(:);

if isfield(Prob,'x_L')
   x_L=Prob.x_L(:);
else
   x_L=[];
end
if isfield(Prob,'x_U')
   x_U=Prob.x_U(:);
else
   x_U=[];
end

if isempty(x_L),      x_L =-Inf*ones(n,1); end
if isempty(x_U),      x_U = Inf*ones(n,1); end
if isempty(b_L),      b_L =-Inf*ones(m,1); end
if isempty(b_U),      b_U = Inf*ones(m,1); end

% Necessary?
x_L = [x_L;-Inf*ones(n-length(x_L),1)]; % Change x_L to full length
x_U = [x_U; Inf*ones(n-length(x_U),1)]; % Change x_U to full length
b_L = [b_L;-Inf*ones(m-length(b_L),1)]; % Change b_L to full length
b_U = [b_U; Inf*ones(m-length(b_U),1)]; % Change b_U to full length

% Find which variables and constraints are equalities

xEqual=x_L==x_U;
xFree=find(~xEqual);

bEqual=b_L==b_U;

if TransfType==1

   x=x_L;

   % Change all -Inf lower bounds and very negative bounds to LowInf
   x(x_L <=LowInf)=LowInf;  % Transform so x-x_L >= 0
   ixU=x_U <=abs(LowInf) & ~xEqual;      % Check which upper bounds on x to use

   % Check which lower bounds on lhs that should be used
   %ibL=~isinf(b_L) | b_L >= -1E6 & ~bEqual;            
   ibL=b_L >= -1E6 & ~bEqual;            

   % Check which upper bounds on rhs that should be used
   %ibU=~isinf(b_U) | b_U <=  1E6 & ~bEqual;            
   ibU= b_U <=  1E6 & ~bEqual;            

   me=sum(bEqual);
   mL=sum(ibL);
   mU=sum(ibU);
   mxU=sum(ixU);
   nF=length(xFree);
   Ax=A*x;

   bL=b_L(find(ibL))-Ax(find(ibL));
   bU=b_U(find(ibU))-Ax(find(ibU));

   iPosL=bL > 0; % These must be put as equalities
   iPosU=bU < 0; % These must be put as equalities

   mPL=sum(iPosL);
   mPU=sum(iPosU);
   meq=me+mPL+mPU;
   AA=zeros(me+mL+mU+mxU,nF+mPL+mPU);
   mAA=size(AA,1);
   bb=zeros(me+mL+mU+mxU,1);

   % First me equality constraints. Must have positive rhs.
   iE=find(bEqual);
   b1=b_L(iE)-Ax(iE);
   iPos1=b1 > 0; % These keep their signs
   nP=sum(iPos1);
   bb(1:me)        = abs(b1);

   AA(1:nP,   1:nF)=  A (iE(find( iPos1)),xFree);
   AA(nP+1:me,1:nF)= -A (iE(find(~iPos1)),xFree);
   bb(1:nP)        =  b1(   find( iPos1));
   bb(nP+1:me)     = -b1(   find(~iPos1));

   % Lower bounds with positive rhs must be equalities, 
   % because slacks has negative coefficients
   % Same thing with Upper bounds and negative rhs

   bb(me+1        :me+mPL)          =  bL(find( iPosL));
   bb(me+mPL+1    :meq)             = -bU(find( iPosU));
   bb(meq+1       :meq+mL-mPL)      = -bL(find(~iPosL));
   bb(meq+mL-mPL+1:me+mL+mU)        =  bU(find(~iPosU));
   bb(me+mL+mU+1  :mAA)             = x_U(find( ixU))-x(find(ixU));

   zL=find(ibL);
   zU=find(ibU);
   AA(me+1        :me+mPL    ,1:nF) =  A (zL(find( iPosL)),1:nF);
   AA(me+mPL+1    :meq       ,1:nF) = -A (zU(find( iPosU)),1:nF);
   AA(meq+1       :meq+mL-mPL,1:nF) = -A (zL(find(~iPosL)),1:nF);
   AA(meq+mL-mPL+1:me+mL+mU  ,1:nF) =  A (zU(find(~iPosU)),1:nF);
   B=eye(nF);
   AA(me+mL+mU+1  :mAA       ,1:nF) =  B (find(ixU),:);

   % Set part corresponding to slack variables with negative coefficients,
   % surplus variables.

   AA(me+1:me+mPL+mPU,nF+1:nF+mPL+mPU)=-eye(mPL+mPU);

   % New simplex problem
   % min c'(x-x_L) s/t
   %
   %  A (x-x_L)             = b_L - A x_L for equalities
   %                                     Inequalities:
   %  A (x-x_L) - I s_1     = bL = b_L - A x_L for bL > 0 . s_1 slack vars
   % -A (x-x_L) - I s_2     =-bU = b_U - A x_L for bU < 0 . s_2 slack vars
   % -A (x-x_L)            <=-bL = b_L - A x_L for bL <=0 
   %  A (x-x_L)            <= bU = b_U - A x_L for bU >=0 
   %  I (x-x_L)            <= x_U - x_L        for upper bounds < 1E4 (Inf)
   %                                           and x_U > x_L (free variables)
   % x-x_L,s_1,s_2 >= 0

   % Translation:
   %   if nF~=0  % No free variables
   %      x(xFree)=y(1:nF)+x(xFree); % New x = x + x_L
   %   end

elseif TransfType==2 
   eqCon   = find( bEqual & ~isinf(b_L));
   ineqLow = find(~bEqual & ~isinf(b_L));
   ineqUpp = find(~bEqual & ~isinf(b_U));
   % Put linear equalities first in A and 
   % turn '<=' inequalities to '>=' inequalities
      
   AA  = [A(eqCon,:);-A(ineqLow,:);A(ineqUpp,:)];
   bb  = [b_L(eqCon);-b_L(ineqLow);b_U(ineqUpp)];
   meq = length(eqCon);

elseif TransfType==3
   eqCon   = find( bEqual & ~isinf(b_L));
   ineqLow = find(~bEqual & ~isinf(b_L));
   ineqUpp = find(~bEqual & ~isinf(b_U));

   % Put linear equalities first in A, followed by fixed variables and 
   % turn '<=' inequalities to '>=' inequalities
   % Add upper bounds last
   iFix = find(~isinf(x_U) &  xEqual);
   iUpp = find(~isinf(x_U) & ~xEqual);
   I    = eye(n,n);
      
   AA  = [A(eqCon,:);I(iFix,:);-A(ineqLow,:);A(ineqUpp,:);I(iUpp,:)];
   bb  = [b_L(eqCon);x_U(iFix);-b_L(ineqLow);b_U(ineqUpp);x_U(iUpp)];
   meq = length(eqCon) + length(iFix);
end

if makeEQ
   m=size(AA,1);
   mI = m-meq;
   AA=[AA, [zeros(meq,mI);eye(mI,mI)]];
end

% MODIFICATION LOG
%
% 990617  hkh  Modifications for MIDEVA
% 090813  med  mlint check