% Prints vector v in the format
%
% Name("1st item",    "item vRow")        = v_1 v_2 ... v_vRow
% Name("item vRow+1", "item vRow+2*vRow") = v_vRow+1 v_vRow+2 ... v_2*vRow
% etc.
%
% where the actual numbers in " " are computed and displayed. 
%
% function vPrint(v,Name,FORM,vRow)
%
% v    Row or column vector. For v matrix, v=v(:) is used.
% Name Name of vector. Default 'x'.
% FORM Optional format string. Default ' %7.2f'.
% vRow Number of elements per row

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written June 18, 1999. Last modified June 18, 1999.

function vPrint(v,Name,FORM,vRow)

if nargin < 4
   vRow=[];
   if nargin < 3
      FORM=[];
      if nargin < 2
         Name='x';
      end
   end
end

if isempty(FORM), FORM=' %7.2f'; end

v=v(:);
Nv=length(v);

s=sprintf('(%d:%d)',Nv,Nv);
Ns=length(s);

if isempty(vRow)
   global MAXCOLS
   if isempty(MAXCOLS),MAXCOLS=80; end
   vRow=floor((MAXCOLS-(1+Ns+length(Name)))/8);
end

space=' ';
i=1;
while i <= Nv
   if rem(i,vRow)==1 | vRow==1
      fprintf('\n');
      fprintf(Name); 
      s=sprintf('(%d:%d)',i,min(i+vRow-1,Nv));
      fprintf('%s',s); 
      if length(s) < Ns 
         fprintf('%s',space(ones(1,Ns-length(s))));
      end
      fprintf(' '); 
   end
   fprintf(FORM,v(i)); 
   i=i+1;
end
fprintf('\n');