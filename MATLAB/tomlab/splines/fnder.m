%FNDER Differentiate a function.
%
%   FNDER(F)
%
%   returns the (representation of the) first derivative of the
%   univariate function contained in F (and in the same form).  
%
%   FNDER(F,DORDER) returns the DORDER-th derivative, with DORDER expected
%   to be of the form [d1,...,dm] in case the function in F is m-variate,
%   and, for each i=1,..,m,  di  an integer to indicate that the function
%   in F is to be differentiated di-fold with respect to its i-th argument.
%   Here, di may be negative, resulting in di-fold integration with respect
%   to the i-th argument.
%
%   For example,
%
%   fnval( fnder( sp, 2), 3.14 );
% 
%   gives the value at 3.14 of the function in sp, while
%
%   sp0 = fnint( fnder( sp ) );
%
%   gives a function which differs from sp by some constant only (namely, by
%   its value at 0).
%
%   See also FNDIR, FNINT.

function fprime = fnder(f,dorder)


if nargin <2
   dorder=[];
   if nargin <1
      f=[];
   end
end

disp('Dummy routine fnder for TOMLAB');

fprime=0;

