% Run menu program asking the question Q, with the menu options as
% rows in the string matrix S
%
% function [n] = strmenu (Q,S);
%
% If the number of items > MAXMENU (global var, default 50), then
% strmenu asks for the part to choose items from.
%
% The MATLAB menu routine does not take matrices as input,
% therefore this routine makes a long calling string and calls menu
% using the eval function
%
% Input:
%         Q   The menu question as a string
%         S   The different options
%
% OUTPUT: 
%         n        Choice by user
%

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Feb 11, 1998.  Last modified June 21, 1999.

function [n] = strmenu(Q,S)

global MAXMENU
if isempty(MAXMENU)
   MAXMENU=50;
end

m=size(S,1);

if m == 1
   n=1;  % Nothing to choose between
elseif m <= MAXMENU
   %T=['n=menu(''' deblank(Q) ''''];
   %for i=1:size(S,1)
   %    T=[T ',''' deblank(S(i,:)) ''''];
   %end
   %T=[T ');'];
   %eval(T);

   T=cell(size(S,1),1);
   for i=1:size(S,1)
       T{i,1}=S(i,:);
   end
   n=menu(deblank(Q),T);
else

   p=ceil(m/MAXMENU);

   Q2=['. Items > ' num2str(MAXMENU) '. Select part:'];
   S2=cell(p,1);
   k=1;
   for i=1:p-1
       S2{i}=[deblank(S(k,:)) ' - ' deblank(S(k+MAXMENU-1,:))];
       k=k+MAXMENU;
   end
   S2{p}=[deblank(S(k,:)) ' - ' deblank(S(m,:))];
   n=[];
   while isempty(n)
      np=menu([deblank(Q) Q2],S2);

      i0=(np-1)*MAXMENU+1;
      i1=min(i0+MAXMENU-1,m);
      T=cell(i1-i0+2,1);
      for i=i0:i1
          T{i-i0+1}=deblank(S(i,:));
      end
      T{i1-i0+2}=['Select another part of menu items '];
      n=menu(deblank(Q),T);
      if n==i1-i0+2, n=[]; end
   end
   n=(np-1)*MAXMENU+n;
end