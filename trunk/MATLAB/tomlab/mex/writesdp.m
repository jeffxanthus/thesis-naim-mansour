%  This function takes a problem in SeDuMi MATLAB format and writes it out 
%  in SDPpack format.  
%
%  Usage:
%
%  writesdp(fname,A,b,c,K)
%
%      fname           Name of SDPpack file, in quotes
%      A,b,c,K         Problem in SeDuMi form
%
%  Notes:
%
%     Problems with complex data are not allowed.    
%
%     Rotated cone constraints are not supported.  
%
%     Nonsymmetric A.s and C.s matrices are symmetrized with A=(A+A')/2
%     a warning is given when this happens.
%
%     Floating point numbers are written out with 18 decimal digits for
%     accuracy.
%
%  Please contact the author (Brian Borchers, borchers@nmt.edu) with any
%  questions or bug reports.
% 
function writesdp(fname,A,b,c,K)
%
%  First, check for complex numbers in A, b, or c.
%
if (isreal(A) ~= 1),
  fprintf('A is not real!\n');
  return
end
if (isreal(b) ~= 1),
  fprintf('b is not real!\n');
  return
end
if (isreal(c) ~= 1),
  fprintf('c is not real!\n');
  return
end
%
%  Check for any rotated cone constraints.
%
if (isfield(K,'r') & (~isempty(K.r)) & (K.r ~= 0)),
  fprintf('rotated cone constraints are not yet supported.\n');
  return
end
%
%  Get the size data.
%
if (isfield(K,'l')),
  nlin=K.l;
  sizelin=nlin;
  if (isempty(sizelin)),
    sizelin=0;
    nlin=0;
  end
  if (K.l == 0),
    nlin=0;
    sizelin=0;
  end
else
  nlin=0;
  sizelin=0;
end

if (isfield(K,'s')),
  nsdpblocks=length(K.s);
  sizesdp=sum((K.s).^2);
  if (isempty(sizesdp)),
    sizesdp=0;
    nsdpblocks=0;
  end
  if (K.s == 0),
    nsdpblocks=0;
    sizesdp=0;
  end
else
  sizesdp=0;
  nsdpblocks=0;
end

if (isfield(K,'q')),
  nqblocks=length(K.q);
  sizeq=sum(K.q);
  if (isempty(sizeq)),
    sizeq=0;
    nqblocks=0;
  end
  if (K.q == 0),
    nqblocks=0;
    sizeq=0;
  end
else
  nqblocks=0;
  sizeq=0;
end
%
%  Find the number of constraints.
%
m=length(b);
%
%  Deal with the following special case.  If A is transposed, transpose
%  it again so that it is of the right size.
%
[Am,An]=size(A);
if (Am ~= m),
  if (An == m),
    fprintf('Transposing A to match b \n');
    A=A';
  else
    fprintf('A is not of the correct size to match b \n');
    return;
  end
end
%
%  Deal with the following special case:  if c==0, then c should really
%  be a zero vector of the appropriate size.
%
if (c == 0),
  fprintf('Expanding c to the appropriate size\n');
  [Am,An]=size(A);
  c=zeros(1,An);
end
%
% SeDuMi uses min problems and SDPpack uses max, so we need to flip the
% sign of the objective function. 
%
c=-c;
%
%  print out some size information
%
fprintf('Number of constraints: %d \n',m);
fprintf('Number of SDP blocks: %d \n',nsdpblocks);
fprintf('Number of QC blocks: %d \n',nqblocks);
fprintf('Number of LP vars: %d \n',nlin);
%
%  Open up the file for writing.
%
fid=fopen(fname,'w');
%
%  Print out m, the number of constraints.
%
fprintf(fid,'%d \n',m);
%
%  Next, b, with one entry per line.
%
fprintf(fid,'%.18e\n',full(b));
%
%  Next, the semidefinite part.
%
if (nsdpblocks == 0),
  fprintf(fid,'0\n');
else
%
%  Print out the number of semidefinite blocks.
%
  fprintf(fid,'%d\n',nsdpblocks);
%
%  For each block, print out its size.
%
  fprintf(fid,'%d\n',full(K.s));
%
%  Next, the cost matrix C.s.
%
%
%  First, calculate where in c things start.
%
  base=sizelin+sizeq+1;
%
%  Next, work through the blocks.
%
  for i=1:nsdpblocks,
    fprintf(fid,'1\n');
    work=c(base:base+K.s(i)^2-1);
    work=reshape(work,K.s(i),K.s(i));
    if (work ~= work')
      fprintf('Non symmetric C.s matrix!\n');
      work=(work+work')/2;
    end
    work=triu(work);
    [II,JJ,V]=find(work);
    cnt=length(II);
    fprintf(fid,'%d\n',cnt);
    if (cnt ~= 0)
      fprintf(fid,'%d\n%d\n%.18e\n',[II JJ V]');
    end
%
%  Next, update to the next base.
%
    base=base+K.s(i)^2;
  end
%
%  Now, loop through the constraints, one at a time.
%
  for cn=1:m,
%
%  Print out the SDP part of constraint cn.
%
    base=sizelin+sizeq+1;
    for i=1:nsdpblocks,
      fprintf(fid,'1\n');
      work=A(cn,base:base+K.s(i)^2-1);
      work=reshape(work,K.s(i),K.s(i));
      if (work ~= work'),
        fprintf('Non symmetric A.s matrix! \n');
        work=(work+work')/2;
      end
      work=triu(work);

      [II,JJ,V]=find(work);
      cnt=length(II);
      fprintf(fid,'%d\n',cnt);
      if (cnt ~= 0),
        fprintf(fid,'%d\n%d\n%.18e\n',[II JJ V]');
      end
%
%  Next, update to the next base.
%
      base=base+K.s(i)^2;
    end

% Done with constraint cn
%
  end
%
% Done with SDP part.
%
end
%
% Next, handle the Quadratic part.
%
%
% Describe the Q blocks.
%
if (nqblocks == 0),
  fprintf(fid,'0\n');
else
  fprintf(fid,'%d\n',nqblocks);
  fprintf(fid,'%d\n',full(K.q));
%
%  Find C.q.
%
  base=sizelin+1;
  cq=c(base:base+sizeq-1);
%
% Print out the C.q coefficients.
%
  fprintf(fid,'%.18e\n',full(cq));
%
%  Next, the constraint matrix A.q. 
%
  Aq=A(:,base:base+sizeq-1);
%
%  Print out the count of nonzeros.
%
  [II,JJ,V]=find(Aq);
  cnt=length(II);
  fprintf(fid,'1\n');
  fprintf(fid,'%d\n',cnt);
  if (cnt ~= 0),
    fprintf(fid,'%d\n%d\n%.18e\n',[II JJ V]');
  end
%
% End of handling quadratic part.
%
end
%
%
% Finally, handle the linear part.
%
if (nlin == 0),
  fprintf(fid,'0\n');
else
%
% Print out the number of linear variables.
%
  fprintf(fid,'%d\n',nlin);
%
% Print out C.l
%
  fprintf(fid,'%.18e\n',full(c(1:nlin)));
%
%  Print out the A matrix.
%
  Al=A(:,1:nlin);
  [II,JJ,V]=find(Al);
  cnt=length(II);
  fprintf(fid,'1\n');
  fprintf(fid,'%d\n',cnt);
  if (cnt ~= 0),
    fprintf(fid,'%d\n%d\n%.18e\n',[II JJ V]');
  end
end