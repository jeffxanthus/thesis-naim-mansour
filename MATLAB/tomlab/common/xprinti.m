% xprinti.m: 
% Print the real vector x in integer format, 
% using the number of digits given by the variable Dig.
% If Dig is empty, the number of digits are determined by xprinti.
%
% nRow elements are printed on each row. If empty, the maximal
% number possible is determined by xprinti.
% xprinti assumes 80 columns if not the global variable MAXCOLS is set.
%
% On the first row, the Name is written (if not empty)
% The other rows have length(Name) spaces first
%
% function xprinti(x,Name,Dig,nRow)
%
% x    Row or column vector. For x matrix, x=x(:) is used.
% Name Name of vector. Default ''.
% Dig  Number of digits. Format ' %Dig.0f' used.
% nRow Number of elements per row. 

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 1997-2005 by Tomlab Optimization Inc., $Release: 4.7.0 $
% Written Aug 15, 1997.   Last modified Jan 17, 2005.

function xprinti(x,Name,Dig,nRow)

global MAXCOLS

if nargin < 4
   nRow=[];
   if nargin < 3
      Dig=[];
      if nargin < 2
         Name=[];
      end
   end
end

if isempty(Dig)
   x=full(x(:));
   maxx=max(x);
   minx=min(x(find(x<0)));
   if isempty(minx), minx=-1; end
   if maxx > 0
      Dig=real(max(ceil(log10(maxx)),1+ceil(log10(minx))));
   else
      Dig=real(1+ceil(log10(minx)));
   end
end

FORM=sprintf(' %%%d.0f',Dig) ;

if isempty(nRow)
   if isempty(MAXCOLS),MAXCOLS=80; end
   nRow=floor((MAXCOLS-(1+length(Name)))/(Dig+1));
end

xprint(x,Name,FORM,nRow);

% MODIFICATION LOG:
%
% 970815 hkh  Written
% 050117 med  mlint revision