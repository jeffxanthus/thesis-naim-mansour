% Print the real vector x, using the format statement FORM.
% nRow elements are printed on each row.
% On the first row, the Name is written (if not empty)
% The other rows have length(Name) spaces first
%
% function xprint(x,Name,FORM,nRow)
%
% x    Row or column vector. For x matrix, x=x(:) is used.
% Name Name of vector. Default ''.
% FORM Optional format string. Default ' %10.6f'.
% nRow Number of elements per row

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 1995-2005 by Tomlab Optimization Inc., $Release: 4.7.0 $
% Written Mar 1, 1995.    Last modified Jan 17, 2005.

function xprint(x,Name,FORM,nRow)

global MAXCOLS

if nargin < 4
   nRow=6;
   if nargin < 3
      FORM=[];
      if nargin < 2
         Name=[];
      end
   end
end

if isempty(FORM), FORM=' %10.6f'; end

x=full(x(:));
n=length(x);
N=length(Name);
if isempty(nRow)
   if isempty(MAXCOLS),MAXCOLS=80; end
   nRow=floor((MAXCOLS-(1+N))/11);
end

if ~isempty(Name), fprintf(Name); end
space=' ';
i=1;
while i <= n
   if rem(i,nRow)==1 & i > 1
      fprintf('\n');
      if ~isempty(Name), fprintf('%s',space(ones(1,N))); end
   end
   fprintf(FORM,x(i)); 
   i=i+1;
end
fprintf('\n');

% MODIFICATION LOG:
%
% 950301 hkh  Written
% 050117 med  mlint revision