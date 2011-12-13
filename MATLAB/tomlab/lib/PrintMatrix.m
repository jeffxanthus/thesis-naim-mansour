% PrintMatrix prints matrix with labels and tries to be intelligent
%
% function PrintMatrix(A,Name,MaxCols,ColScr,RowLab,ColLab)
%
% A       Matrix
% Name    Name of matrix A
% MaxCols Maximal number of columns per row displayed. If empty, using
%         global MAXCOLS, if defined.
% ColScr  Number of matrix columns per row displayed
% RowLab  Row Label. String var. Separate each label with space.
% ColLab  Column Label. String var. Separate each label with space.
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Aug 11, 1997.   Last modified Aug 13, 2009.

function PrintMatrix(A,Name,MaxCols,ColScr,RowLab,ColLab)

global MAXCOLS
[m,nc] = size(A);

if m==0 || nc==0   % Empty matrix case
   % Empty matrix
   fprintf('PrintMatrix: Empty Matrix ');
   if exist('Name','var'),fprintf('%s',Name); end 
   fprintf('\n\n');
   return
end

space = ' ';

if nargin < 6, 
   ColLab=[]; 
   if nargin < 5
      RowLab=[]; 
      if nargin < 4
         ColScr=[]; 
         if nargin < 3
            MaxCols=[]; 
            if nargin < 2
               Name=[]; 
               if nargin < 1
                  A=[]; 
               end
            end
         end
      end
   end
end

if isempty(ColScr), ColScr=16; end
if isempty(MaxCols)
   if isempty(MAXCOLS)
      MaxCols=80;
   else
      MaxCols=MAXCOLS;
   end
end

if isempty(RowLab)
   RowLab = []; 
   for i=1:min(99,m), RowLab = [RowLab, '--',int2str(i),'--> ']; end
   for i=100:m, RowLab = [RowLab, '-',int2str(i),'-->  ']; end
   RowLab=RowLab(1:length(RowLab)-1);
end

if isempty(ColLab)
   ColLab = [];
   for i=1:nc, ColLab = [ColLab, '----',int2str(i),'---- ']; end
   ColLab=ColLab(1:length(ColLab)-1);
end

RLen=7;  % Length of Row names

z=A(:);
maxA=max(z);
minA=min(z(find(z<0)));
if isempty(minA), minA=-1; end
minA0=min(abs(z(find(z~=0))));
if isempty(minA0), minA0=0; end

if maxA > 0
   l=real(max(ceil(log10(maxA)),1+ceil(log10(minA))));
elseif maxA==0 & minA==0
   l=1;
else
   l=real(1+ceil(log10(minA)));
end
if minA0 <=0
   d=1;
else
   d=full(real(min(abs(floor(log10(minA0))),5)));
end
if d==0
   % Check if some elements need decimals
   if ~all(A==round(A*1000)/1000)
      d=3;
   elseif ~all(A==round(A*100)/100)
      d=2;
   elseif ~all(A==round(A*10)/10)
      d=1;
   end
end

Len=max(6,l+d+2);
if Len > 12 % Allow max length of 12 chars for each item
   l=Len-2-d;
   Len=12;
   MAXA=10^l;
   MINA=-10^(l-1);
else
   MAXA=Inf;
   MINA=-Inf;
   if minA0 < 0.0001
      % Must use exp format, At least length 9 for exp format
      Len=max(9,Len);
   end
end
FORM1=sprintf(' %%%d.2e',Len);     % Use only 2 decimals in exp format
FORM2=sprintf(' %%%d.%df',Len,d);

% Check how many items is possible to display in each row
items=floor((MaxCols-RLen-1)/(Len+1));
if items < ColScr, ColScr=items; end

% Remove extra spaces (delimiters)
ndx1 = find(ColLab==' ');
ndx2 = find([ndx1,0]==[-1,ndx1+1]);
if ~isempty(ColLab), ColLab(ndx1(ndx2))=[]; end

ndx1 = find(RowLab==' ');
ndx2 = find([ndx1,0]==[-1,ndx1+1]);
if ~isempty(RowLab), RowLab(ndx1(ndx2))=[]; end

% Determine position of delimiters

cpos = find(ColLab==' ');
if length(cpos)<nc-1, error('Not enough column labels.'); end
cpos = [0,cpos,length(ColLab)+1];

rpos = find(RowLab==' ');
if length(rpos)<m-1, error('Not enough row labels.'); end
rpos = [0,rpos,length(RowLab)+1];

col=1;
n = min(ColScr-1,nc-1);
%disp(' ')
if ~isempty(Name)
   fprintf('%s',Name);  % Write name on separate row
   fprintf('\n'); 
   s = space(ones(1,RLen+2));
   %if length(Name)+2 < RLen 
   %   s = [Name space(ones(1,RLen+1-(length(Name)+2)))];
   %   %fprintf('%s',space(ones(1,RLen-length(Name)+2)));
   %else
   %   fprintf('\n'); % If Name too long, write on separate row
   %   s = space(ones(1,RLen+2));
   %   fprintf('%s',Name);
   %   fprintf('\n');
   %end
else
   s = space(ones(1,RLen+2));
end

while col<=nc
  if col~=1
     s = space(ones(1,RLen+2));
  end
  % Print labels
  for j=0:n,
      lab = ColLab(cpos(col+j)+1:cpos(col+j+1)-1);
      if length(lab) > Len,
        lab=lab(1:Len);
      else
        lab=[space(ones(1,Len-length(lab))),lab]; 
      end
      s= [s,' ',lab];
  end
  %disp(setstr(s))
  for i=1:m,
      s = RowLab(rpos(i)+1:rpos(i+1)-1); 
      if length(s) > RLen 
         s=s(1:RLen); 
      else 
         s=[space(ones(1,RLen-length(s))),s]; 
      end
      s = [' ',s];
      for j=0:n
          element = A(i,col+j);
          if element==0,
             s=[s,[space(ones(1,Len)),'0']];
          %elseif (element>=1.e6)|(element<=-1.e5)|(abs(element)<.0001)
          elseif (element>=MAXA) | (element<=MINA) | (abs(element)<.0001)
             %s=[s,sprintf(' %12.5e',element)];
             s=[s,sprintf(FORM1,element)];
          else
             %s=[s,sprintf(' %12.5f',element)];
             s=[s,sprintf(FORM2,element)];
          end
      end
      disp(s)
      %fprintf('length s = %d\n',length(s))
  end % for
  col = col+ColScr;
  disp(' ')
  if (nc-col<n), n=nc-col; end
end % while

% MODIFICATION LOG
%
% 981108  hkh  Add checks if d==0, to get up to 3 decimals
% 031216  hkh  Use full to avoid sparse 1x1 d matrix
% 090813  med  mlint check