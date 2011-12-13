%FNVAL Evaluate a function.
%
%   FNVAL(F,X)  evaluates the function  f  in F at the specified points. 
%
%   The output is a matrix of size [d*m,n] if the function in F is
%   univariate, d-vector valued and [m,n] = size(X) .
%
%   If the function in F is d-vector valued and m-variate with m>1, then
%
%                         [d*m,n],       if X is of size [m,n]
%   the output is of size [d,n1,...,nm], if d>1  and X is {X1,...,Xm} 
%                         [n1,...,nm],   if d==1 and X is {X1,...,Xm} 
%
%
%   See also SPVAL, PPUAL, PPVAL.


function values = fnval(f,x)


if nargin <2
   x=[];
   if nargin <1
      f=[];
   end
end

disp('Dummy routine fnval for TOMLAB');

values=0;

