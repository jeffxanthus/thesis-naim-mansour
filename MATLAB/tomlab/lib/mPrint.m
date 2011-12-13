% Prints matrix A in the format
% Name(i,:) a_i1 a_i2 ... a_in
%
% function mPrint(A,Name,FORM,Arow)
%
% A      Matrix
% Name   Name of matrix. Default A.
% FORM   Optional format string. Default ' %7.2f'.
% Arow   Number of elements per row

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Aug 11, 1997.  Last modified June 18, 1999.

function mPrint(A,Name,FORM,Arow)

if nargin < 4
   Arow=[];
   if nargin < 3
      FORM=[];
      if nargin < 2
         Name='A';
      end
   end
end

global MAXCOLS

if isempty(FORM), FORM=' %7.2f'; end

[m,n]=size(A);
s=sprintf('(%0.0f,:)',m);
Ns=length(s);

if isempty(Arow)   
   if isempty(MAXCOLS),MAXCOLS=80; end
   Arow=floor((MAXCOLS-(1+Ns+length(Name)))/8);
end

space=' ';
for j=1:m
    i=1;
    while i <= n
       if rem(i,Arow)==1 | Arow==1
          if i==1
             fprintf('\n');
             fprintf('%s',Name); 
             s=sprintf('(%0.0f,:)',j); 
             fprintf('%s',s); 
             if length(s) < Ns, fprintf('%s',space(ones(1,Ns-length(s)))); end
          else
             fprintf('\n');
             fprintf(Name); 
             fprintf('%s',space(ones(1,Ns)));
          end
          fprintf(' '); 
       end
       fprintf(FORM,A(j,i)); 
       i=i+1;
    end
end
fprintf('\n');