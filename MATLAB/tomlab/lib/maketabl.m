% function maketabl(A, frame, Tsize, lrc, finame, decs, caption, head, ...
%                   subhead, lStr, rStr, Note, hcols, tPlace, rowdecs)
%
% maketabl. Make tex file with floating table environment from matrix A
%           Label is set the same as the name of the file 
%           The table is centered on the page (and floating)
%
% A       Matrix, m by n. If value NaN, nothing is displayed
%         If A is structure, A.text is displayed.
% frame   Type of table.
%         0 = No frame
%         1 = No frame, lines above/below the header and at the end
%         2 = No frame, lines above/below the header and above/below last line
%         3 = Frame around the table
%         4 = Also frame around the header
%         5 = All entries are framed
% Tsize   Type size of table. Tsize < 0 ==> Set no size, use LaTeX default.
%         0 = tiny              1 = scriptsize (default)    2 = footnotesize
%         3 = small             4 = normalsize              5 = large 
% lrc     Adjustment for each column element. Vector is expanded to correct
%         length
%         = 1 left  adjusted
%         = 2 right adjusted
%         = 3 centered
% finame  Name of tex file (without extension .tex)
% decs    Output col j in A with decs(j) decimals. decs(j) < 0 gives integer.
%         Vector is expanded to correct length. Default is integer output.
%         Also see parameter rowdecs 
% caption Table header text
% head    Header. If empty, no header is displayed   
% subhead Cell array with sub-header strings. size(subhead,2)==size(A,2) 
% lStr    If not empty, use this string matrix as first column, to the left
%         size(lStr,1)==size(subhead,1) + size(A,1))
% rStr    If not empty, use this string matrix as last column, to the right
%         size(rStr,1)==size(subhead,1) + size(A,1))
% Note    Three dimensional matrix (m,n,# of notes). If Note(i,j,k)=1,then
%         the number k is set as a raised note (in \tiny) for element(i,j)
% hcols   Number of columns to use. If hcols has length 2, Use
%         columns hcols(1)-hcols(2) for header. NOT IMPLEMENTED YET
% tPlace  Placement of table. h=here,t=top,b=bottom,p=float. Default htbp
% rowdecs 1 iff the figures in decs refers to rows instead of columns. Default 0

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written May 1, 1997.     Last modified Sep 27, 1999.

function maketabl(A, frame, Tsize, lrc, finame, decs, caption, head, ...
                  subhead, lStr, rStr, Note, hcols, tPlace, rowdecs)
               
if nargin < 15 % TH
   rowdecs=0; % TH
   if nargin < 14
      tPlace='htbp';
      if nargin < 13
         hcols=[];
         if nargin < 12
            Note=[];
            if nargin < 11
               rStr=[];
               if nargin < 10
                  lStr=[];
                  if nargin < 9
                     subhead=[];
                     if nargin < 8
                        head=[];
                        if nargin < 7
                           caption=[];
                           if nargin < 6
                              decs=[];
                              if nargin < 5
                                 finame=[];
                                 if nargin < 4
                                    lrc=[];
                                    if nargin < 3
                                       Tsize=[];
                                       if nargin < 2
                                          frame=[];
                                       end
                                    end
                                 end
                              end
                           end
                        end
                     end
                  end
               end
            end
         end
      end
   end
end
               
lrcdef='lrc';
[m,n]=size(A);
if isempty(finame), finame='amat'; end
if isempty(lrc),    lrc=2*ones(n,1); end
if isempty(decs),   decs=-1*ones(n,1); end
if isempty(frame),  frame=5; end
if isempty(Tsize),  Tsize=1; end

TSIZE=str2mat('tiny','scriptsize','footnotesize','small','normalsize','large'); 

lrc=lrc(:);
k=length(lrc);
last=lrc(k);
if length(lrc) < n,lrc=[lrc(1:k);last*ones(n-k,1)]; end

% TH mod 8 lines
k=length(decs);
last=decs(k);
if rowdecs
   if length(decs) < n ,decs=[decs(1:k);last*ones(n-k,1)]; end
else
   decs=decs(:);
   if length(decs) < n,decs=[decs(1:k) last*ones(1,n-k)]; end
end

left=~isempty(lStr);
right=~isempty(rStr);

Fname=[deblank(finame) '.tex'];
%Fname
%pause
f=fopen(Fname,'w');

fprintf(f,'\\begin{table}[%s',tPlace);
fprintf(f,']\n');
if ~isempty(caption)
   fprintf(f,'\\caption{%s',caption);
   fprintf(f,'}\n',caption);
end
fprintf(f,'\\label{%s',deblank(finame));
fprintf(f,'}\n',deblank(finame));
if Tsize >= 0
   fprintf(f,'{\\%s',deblank(TSIZE(min(Tsize+1,size(TSIZE,1)),:)));
   fprintf(f,'\n');
end
if isempty(Note)
   DOUB=0;
else
   DOUB=1; % Just 1pt between 2 columns
   fprintf(f,'\\tabcolsep=1pt');
end
fprintf(f,'\\begin{center}\n');
fprintf(f,'\\begin{tabular}{');
if frame > 2, fprintf(f,'|'); end
if left
   fprintf(f,'l');                  % Left adjusted strings
   if frame > 2, fprintf(f,'|'); end
end
   
for j=1:n
    fprintf(f,'%s',lrcdef(lrc(j)));
    if DOUB, fprintf(f,'@{}l'); end
    if j==n
       if frame > 2, fprintf(f,'|'); end
    else
       if frame > 4, fprintf(f,'|'); end
    end
end
if right
   fprintf(f,'l');                  % Left adjusted strings
   if frame > 2, fprintf(f,'|'); end
end
fprintf(f,'}');
if frame > 0 
   fwrite(f,' \hline'); 
end
fprintf(f,'\n');

if ~isempty(head)
   if left, fprintf(f,'&'); end
   fprintf(f,'\\multicolumn{%d}{',n+DOUB*n);
   if frame > 2, fprintf(f,'|'); end
   fprintf(f,'c');
   if frame > 2, fprintf(f,'|'); end
   fprintf(f,'}\n{{\\bf %s',head);
   fprintf(f,'}}');
   if right, fprintf(f,'&'); end
   fprintf(f,'\n\\\\ ');
   if frame > 0, fwrite(f,' \hline'); end
   fprintf(f,'\n');
end

k=0;
if ~isempty(subhead)
   for i=1:size(subhead,1)
       k=k+1;
       if left
          if i <= size(lStr,1)
             fprintf(f,'%s',deblank(lStr(i,:)));
          end
          fprintf(f,'&\n');
       end
       for j=1:size(subhead,2)
           fprintf(f,'%s',subhead{i,j});
           if DOUB, fprintf(f,'&'); end
           if j<size(subhead,2), fprintf(f,'&'); end
       end
       if right
          fprintf(f,'&');
          if i <= size(rStr,1)
             fprintf(f,'%s',deblank(rStr(i,:)));
          end
       end
       fprintf(f,'\n\\\\ ');
       if frame > 0, fwrite(f,' \hline'); end
       fprintf(f,'\n');
   end
end

for i=1:m
    if left
       if i <= size(lStr,1)
          fprintf(f,'%s',deblank(lStr(k+i,:)));
       end
       fprintf(f,'&\n');
    end
    for j=1:n
        if isstruct(A)
           % Modified to handle structs of strings in A
           fprintf(f,'%s',A(i,j).text); % 
        elseif ~isnan(A(i,j))
           % TH added row/column handling
           if rowdecs
              kk = i;
           else
              kk = j;
           end
           if decs(kk) < 0
              if decs(kk) == -2
                 fprintf(f,'%e',A(i,j)); 
              elseif decs(kk) == -3
                 fprintf(f,'%f',A(i,j));
              elseif decs(kk) == -4
                 fprintf(f,'%g',A(i,j));
              else
                 fprintf(f,'%d',round(A(i,j)));
              end 
           else
              str=['%7.' num2str(decs(kk)) 'f '];
              fprintf(f,str,A(i,j));
           end
        end
        if DOUB
           fprintf(f,'&');
           jj=find(Note(i,j,:));
           % Check if to add any note
           if length(jj) > 0
              fprintf(f,'$^{');
              fprintf(f,'%d',jj(1));
              if length(jj) > 1, fprintf(f,',%d',jj(2:length(jj))); end
              fprintf(f,'}$\n');
           end
        end
        if j~=n, fprintf(f,'&'); end
    end
    if right
       fprintf(f,'&\n');
       if i <= size(rStr,1)
          fprintf(f,'%s',deblank(rStr(k+i,:)));
       end
    end
    fprintf(f,'\n\\\\ ');
    if frame == 2 & i==m-1 
       fwrite(f,' \hline'); 
    end
    if frame > 3 
       fwrite(f,' \hline'); 
    elseif any(frame == (1:3)) & i==m 
       fwrite(f,' \hline'); 
    end
    fprintf(f,'\n');
end

fprintf(f,'\\end{tabular}\n');
fprintf(f,'\\end{center}\n');
if Tsize >= 0
   fprintf(f,'}\n');
end
if DOUB
   fprintf(f,'\\tabcolsep=6pt'); % Reset to normal space between columns
end
fprintf(f,'\\end{table}\n');
fclose(f);

% MODIFICATION LOG
%
% 980321  mbk  A is a matrix of structures containing a string. Used in cls_test
% 980801  th   decs is now either row or column vector,Param 15: rowdecs added
% 981021  hkh  Added JPs handling of negative decimals, combined with THs
%                row/column handling
% 990424  hkh  Add another table format, frame==1, 
% 990916  hkh  Make it work for  Mideva