% -------------------------------------------------------------------------
function [x_L,x_U, Feasible] = preAx(A,b_L,b_U,x_L,x_U,fixV,FV,PriLev)
% -------------------------------------------------------------------------
% -------------------------------------------------
%
% INPUT:
%
% OUTPUT:

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2009 by Tomlab Optimization Inc., $Release: 7.4.0$
% Written Oct 22, 2009.   Last modified Oct 22, 2009.

if nargin < 8
   PriLev = [];
end
[m,n] = size(A);

if isempty(PriLev)
   if n > 10 
      PriLev = 2;
   else
      PriLev = 5;
   end
end

x_L(fixV) = FV;
x_U(fixV) = FV;
PriLev = 2;

AL    = A.*(ones(m,1)*x_L');
AU    = A.*(ones(m,1)*x_U');
Amin  = min(AL,AU);
Amax  = max(AL,AU);
Feasible = 1;

if PriLev > 2
   mPrint(A,'A:   ',[],15)
end
if PriLev > 3
   mPrint(AL,'AL:  ',[],15)
   mPrint(AU,'AU:  ',[],15)
   mPrint(Amin,'Amin:',[],15)
   mPrint(Amax,'Amax:',[],15)
end

Alow  = sum(Amin')';
Aupp  = sum(Amax')';
if PriLev > 1
   disp('[b_L Alow Aupp b_U]');
   disp([b_L full(Alow) full(Aupp) b_U]);
end

% If any of these violated, infeasible choice
Test1 = full(any(Alow > b_U));
Test2 = full(any(Aupp < b_L));
if PriLev > 0
   fprintf('If any of these two violated, infeasible choice:')
   fprintf('Alow > b_U %d ',Test1);
   fprintf('Aupp < b_L %d ',Test2);
   fprintf('\n');
end
if Test1 == 1 | Test2 == 1
   Feasible = 0;
   return
end
if PriLev > 1
   fprintf('\n');
   [x_L x_U]'
end
k = 1;
Change = 1;
while k < 4 & Change == 1
   Change = 0;
   ix = find(Alow == b_U & Aupp ~= b_L);
   if ~isempty(ix)
      if PriLev > 2
         fprintf('Rows were variables could be fixed to lower:')
         xprinti(ix);
      end
      Change = 1;
      for l=1:length(ix)
          i = ix(l);
          j = find(x_L~=x_U & full(A(i,:)' ~= 0));
          if length(j) == 1
             if PriLev > 0
                fprintf('Upper bound on variable %d set to lower bound %f\n',j,x_L(j));
             end
             x_U(j) = x_L(j);
             AU(:,j) = A(:,j)*x_U(j);
          else
             if PriLev > 1
                xprinti(j,'Vars j:');
             end
          end
      end
   end
   iz = find(Aupp == b_L & Alow ~= b_U);
   if ~isempty(iz)
      if PriLev > 2
         fprintf('Rows were variables could be fixed to upper:')
         xprinti(iz);
      end
      Change = 1;
      for l=1:length(iz)
          i = iz(l);
          j = find(x_L~=x_U & A(i,:)' ~= 0);
          if length(j) == 1
             if PriLev > 0
                fprintf('Lower bound on variable %d set to upper bound %f\n',j,x_U(j));
             end
             x_L(j) = x_U(j);
             AL(:,j) = A(:,j)*x_L(j);
          else
             if PriLev > 1
                xprinti(j,'Vars j:');
             end
          end
      end
   end
   if Change
      Amin  = min(AL,AU);
      Amax  = max(AL,AU);
      if PriLev > 2
         mPrint(AL,'AL:  ',[],15)
         mPrint(AU,'AU:  ',[],15)
         mPrint(Amin,'Amin:',[],15)
         mPrint(Amax,'Amax:',[],15)
      end
      Alow  = sum(Amin')';
      Aupp  = sum(Amax')';
      if PriLev > 1
         disp('[b_L Alow Aupp b_U]');
         disp([b_L full(Alow) full(Aupp) b_U]);
         [x_L x_U]'
      end
   end
   k = k+1;
end
if m == 1
   nA = A~=0;
else
   nA = sum(A~=0);
end
if PriLev > 1
   fprintf('\n');
   xprint(nA,'A~0:'),
end
ix = find(nA>0);
for i=ix
    iy = find(A(:,i)~=0);
    v  = A(iy,i);
    bu = (b_U(iy)-Alow(iy)+Amin(iy,i))./v;
    bl = (b_L(iy)-Aupp(iy)-Amax(iy,i))./v;
    if PriLev > 1
       xprint([bu(sign(v)==1);bl(sign(v)==-1)]);
       xprint([bu(sign(v)==-1);bl(sign(v)==1)]);
    end
    %mbu = max([bu(sign(v)==1);bl(sign(v)==-1)]);
    %mbl = min([bu(sign(v)==-1);bl(sign(v)==1)]);
    %xprint(bu,'bu:')
    %xprint(bl,'bl:')
    %fprintf('Var %d: x_L %f ',i,x_L(i));
    %fprintf('mbl %f, mblu %f x_U %f ',mbl,mbu,x_U(i));
    %fprintf('\n');
    mbu = min([bu(sign(v)==1);bl(sign(v)==0)]);
    mbl = max([bu(sign(v)==0);bl(sign(v)==1)]);
    if PriLev > 0
       fprintf('Var %d: x_L %f ',i,x_L(i));
       fprintf('mbl %f, mblu %f x_U %f ',mbl,mbu,x_U(i));
       fprintf('\n');
    end
    if mbu < x_L(i) | mbl > x_U(i)
       Feasible = 0;
    end
    
end

if PriLev > 1
   [x_L x_U]'
end
Feasible

% MODIFICATION LOG
%
% 091022  hkh  Written
