% Converting a general assignment problem (GAP) to general form suitable for
% a MIP solver
%
% Either binary or sos1 variables are used.
%
% function [c, x_L, x_U, b_L, b_U, a, sos1] = abc2gap( A, b, C, SOS1);
%
% GAP problem:
%         m    n
%    min sum  sum c(i,j) * x(i,j)
%        i=1  j=1
%
% subject to
%         n
%        sum x(i,j) = 1   for i=1, ...,m
%        j=1
%
%         m
%        sum a(i,j) * x(i,j) <= b(j)   for j=1, ...,n
%        i=1
%
%        x binary 0/1
%
% INPUT:
% A      m x n constraint matrix
% C      m x n cost matrix
% b      m x 1 right hand side 
% SOS1   If true, generate output for sos1 handling with XPRESS-MP.
%        Otherwise generate output giving an equivalent formulation
%        with standard integer variables
%
% OUTPUT
% c      m * n cost vector
% x_L    m * n lower bounds on decision variables x
% x_U    m * n upper bounds on decision variables x
% b_L    m + n lower bounds on the constraints
% b_U    m + n upper bounds on the constraints
% a      m + n x m * n constraint matrix 
% sos1   Structure with information for xpress.m, MEX-file interface.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., $Release: 4.7.0 $
% Written Oct 27, 1999.   Last modified Jan 17, 2005.

function [c, x_L, x_U, b_L, b_U, a, sos1] = abc2gap( A, b, C, SOS1)

if nargin < 4
   SOS1=0;
end

[m,n]=size(A);

N=m*n;
x_L = zeros(N,1);
x_U = ones(N,1);

a=spalloc(n,N,N);
k=0;
for j=1:n
    a(j,k+[1:m])=A(:,j)';
    k=k+m;
end
b_L = -Inf * ones(n,1);
b_U = b;

c=C(:);

if SOS1
   % Generate the set constraints that defines the sos1 constraints
   A2=zeros(1,N);
   v=m*[0:n-1];
   for j=1:m
       sos1(j).var = j+v; 
       % Define the cost coefficients as the ordering, row 0.
       % Ties are solved by the xpress.m interface.
       sos1(j).row = 0;
       if 0 % NOT use special ordering. 
          sos1(j).row = size(a,1)+1;
          % The ordering info could be c, a, just 1:n
          %A2(1,j+v)=1:n;
          %A2(1,j+v)=A(j,:);
          % Just put the cost coefficients as the ordering.
          % Ties are solved by the xpress.m interface.
          A2(1,j+v)=C(j,:);
       end
       nSum=sum(1:n);
   end
   %if 0 % NOT USE
   %   A2sum=sum(A2);
   %   %b_L = [b_L;-Inf];
   %   %b_U = [b_U;Inf];
   %   b_L = [b_L;min(A2sum,0)];
   %   b_U = [b_U;max(A2sum,0)];
   %   a= [a;A2];
   %end

   % Generate A2, the sos1 constraints
   A2=spalloc(m,N,N);
   v=m*[0:n-1];
   for j=1:m
       A2(j,j+v)=1;
   end
   b_L = [b_L;ones(m,1)];
   b_U = [b_U;ones(m,1)];
   a= [a;A2];
else
   % Generate A2, the sos1 constraints
   A2=spalloc(m,N,N);
   v=m*[0:n-1];
   for j=1:m
       A2(j,j+v)=1;
   end
   b_L = [b_L;ones(m,1)];
   b_U = [b_U;ones(m,1)];
   a= [a;A2];
   sos1=[];
end

% MODIFICATION LOG:
%
% 050117 med  mlint review