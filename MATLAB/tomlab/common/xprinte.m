% Print the real vector x in exponential format, 
% using the format statement FORM.
% nRow elements are printed on each row.
% On the first row, the Name is written (if not empty)
% The other rows have length(Name) spaces first
%
% function xprinte(x,Name,FORM,nRow)
%
% x    Row or column vector. For x matrix, x=x(:) is used.
% Name Name of vector. Default ''.
% FORM Optional format string. Default ' %14.6e'.
% nRow Number of elements per row. 

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 1995-2004 by Tomlab Optimization Inc., $Release: 4.7.0 $
% Written Mar 1, 1995.    Last modified June 18, 1999.

function xprinte(x,Name,FORM,nRow)

if nargin < 4
   nRow=6;
   if nargin < 3
      FORM=' %14.6e';
      if nargin < 2
         Name=[];
      end
   end
end

xprint(x,Name,FORM,nRow);
