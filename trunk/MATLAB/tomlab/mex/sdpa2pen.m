%%*****************************************************************************
% Convert file in SDPA format to PENSDP format
%
% [pen] = sdpa2pen(fname, Dir)
%
% Input: fname...name of the file containing SDP data in SDPA format
% Input: Dir  ...directory name for the SDPA data file (optional) 
%
% Output: pen...structure with data arrays in PENDSP format
%
% Copyright (c) 2002 by M. Kocvara and M. Stingl
% Version 29/4/2002
% Revised by Kenneth Holmstrom July 1, July 3, 2002
%*****************************************************************************

function [pen] = sdpa2pen(fname,Dir)

% Holmstrom revision
if nargin > 1
   Dir = deblank(Dir);
   if ispc
      if Dir(end) == '\'
         fname = [Dir  fname];
      else
         fname = [Dir '\'  fname];
      end
   else
      if Dir(end) == '/'
         fname = [Dir  fname];
      else
         fname = [Dir '/'  fname];
      end
   end
end

fname

if exist(fname)
   fid = fopen(fname,'r');
else
   pen = [];
   fprintf('File %s not found, please enter a valid name.\n',fname);
   return
end

% End of Holmstrom revision

% ignore the comment lines at the file beginning
comments = 0;
while comments == 0
   tmp = fgetl(fid);
   k1=findstr(tmp,'"');   
   k2=findstr(tmp,'*');
   if (isempty(k1) & isempty(k2)),
      comments = 1;
   end
end

% read the problem sizes
vars = str2num(tmp);
pen.vars = vars;

nblc = fgetl(fid);
nbl = str2num(nblc);
iblc = fgetl(fid);
ibl = str2num(iblc);

% read the objective function vector, ignoring special characters
tmp = fgetl(fid);
brac = '{}(),';
tind=[];
for i=1:length(brac);
   ind = findstr(tmp,brac(i));
   tind = [tind,ind];
end;
tmp(tind) = blanks(length(tind));
pen.fobj = str2num(tmp);

% read the problem data
[tmp,count] = fscanf(fid,'%f');
len = length(tmp);
ndum = len/5;
dum = reshape(tmp,5,ndum)';

% convert sdpa problem data into pensdp format
% linear constraints
iicon = 0;
pen.constr = 0; pen.ci = 0; pen.bi_dim = 0; pen.bi_idx = 0; pen.bi_val = 0;
concount = 0;
for i=1:nbl
   if (ibl(i) < 0)
      iicon = iicon + 1;
      pen.constr = -ibl(i);
      pen.ci = zeros(1,pen.constr);
      bi = sparse([],[],[],1,1,0);   
      pen.bi_dim = zeros(1,pen.constr);
      for j=1:ndum
         if(dum(j,2) == i)
            concount = concount + 1;
            if(dum(j,1) == 0)
               pen.ci(dum(j,3)) = -dum(j,5);
            else
               bi(dum(j,3),dum(j,1)) = -dum(j,5);
               pen.bi_dim(dum(j,3)) = pen.bi_dim(dum(j,3)) + 1;
            end
         end
      end
   end
end

if (pen.constr > 0 )
   pen.bi_idx = zeros(1,sum(pen.bi_dim));
   pen.bi_val = zeros(1,sum(pen.bi_dim));
   counter = 1;
   for iib = 1:pen.constr
      [ibi,jbi,vbi] = find(bi(iib,:));
      for(kkk = 1 : length(jbi))
         pen.bi_idx(counter) = jbi(kkk)-1;
         pen.bi_val(counter) = vbi(kkk);
         counter = counter + 1;
	   end   
   end
end

% matrix constraints
pen.mconstr = nbl - iicon;
imcon = 0;
pen.ai_dim = zeros(1,pen.mconstr);
pen.ai_row = zeros(1,ndum-concount); pen.ai_col = zeros(1,ndum-concount); 
pen.ai_val = zeros(1,ndum-concount);
ai_nzs = sparse([],[],[],pen.mconstr,pen.vars+1,0);
ijk = 0;
for i=1:nbl
   if (ibl(i) > 0)
      imcon = imcon + 1;
      pen.msizes(imcon) = ibl(i);  
      for j=1:ndum
         if(dum(j,2) == i)
            ai_nzs(imcon,dum(j,1)+1) = ai_nzs(imcon,dum(j,1)+1) + 1;
            ijk = ijk+1;
            pen.ai_row(ijk) = dum(j,3)-1;
            pen.ai_col(ijk) = dum(j,4)-1;
            if(dum(j,1) == 0)
               pen.ai_val(ijk) = dum(j,5);
            else
               pen.ai_val(ijk) = -dum(j,5);
            end
         end
      end
   end
end

[iai,jai,vai] = find(ai_nzs);
[idum,iso]=sort(iai);
pen.ai_idx = jai(iso)'-1;
pen.ai_nzs = vai(iso)';
for i=1:pen.mconstr
   pen.ai_dim(i) = length(iai)-nnz(iai-i);
end         

pen.x0 = zeros(1,vars);
pen.ioptions(1) = 0;
pen.foptions = [];

fclose(fid);