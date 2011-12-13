%% Generator of random BMI (and quadratic MI) problems of the type
%%
%% min lambda    (lambda = x_n)
%% s.t.
%% x_i \in [a,b] i=1,..,n-1   (random bounds)
%% A_0 + sum_{i=1}^{n-1} x_i A_i + sum_{i=1}^{n-1}sum_{j=1}^{n-1} x_i x_j K_{ij} - Id.lambda =< 0
%%
%% Copyright (c) 2002 by M. Kocvara and M. Stingl, PENOPT GbR
%% Version 01/11/2002
%% 
%% Modified to function mode by Anders Goran, Tomlab Optimization Inc. 
%%
%% function pen = bmiran(n,m,dens)
%%
%% n    - number of variables
%% m    - size of matrices
%% dens - density of matrices

function pen = bmiran(n,m,dens)

if nargin < 3
   dens = 0.1;
   if nargin < 2
      m = 5;
      if nargin < 1
         n = 5;
      end
   end
end

constr = 2*(n-1);

pen.vars = n;
pen.constr = constr;
pen.mconstr = 1;
pen.msizes = m;
pen.x0=zeros(1,n);
pen.fobj = [zeros(1,(n-1)),1];
pen.ci= 100.*abs(randn(1,constr));
pen.bi_dim=ones(1,constr);
for i=1:n-1
   pen.bi_idx(1,2*i-1) = i-1;
   pen.bi_val(1,2*i-1) =-1;
   pen.bi_idx(1,2*i)   = i-1;
   pen.bi_val(1,2*i)   = 1;
end

pen.ai_dim=n+1;
pen.ai_idx=(0:n)';

kk = 1;
for i=1:n
   a = sprandn(m,m,dens);
   a = triu(a);
   [ii,jj,ss]=find(a);
   pen.ai_nzs(i) = length(ss);


   ll = kk + length(ss)-1;
   pen.ai_val(kk:ll) = ss;
   pen.ai_row(kk:ll) = ii-1;
   pen.ai_col(kk:ll) = jj-1;
   kk = kk + length(ss);
end


   pen.ai_nzs(n+1) = m;
   pen.ai_val(kk:kk+m-1) = -ones(1,m);
   pen.ai_row(kk:kk+m-1) = 0:m-1;
   pen.ai_col(kk:kk+m-1) = 0:m-1;

n1 = n - 1;% n1 = 2;
pen.ki_dim = n1*(n1-1)/2;  % pen.ki_dim=[n1*(n1+1)/2];
pen.ki_idx = zeros(n1*(n1-1)/2,1);
pen.kj_idx = zeros(n1*(n1-1)/2,1);
kk = 1;
kcount = 1;
for i=1:n1
   for j=i+1:n1
  % for j=i:n1
      pen.ki_idx(kcount)=i;
      pen.kj_idx(kcount)=j;


      a = sprandn(m,m,dens);
      a = triu(a);
      [ii,jj,ss]=find(a);
      pen.ki_nzs(kcount) = length(ss);
      kcount = kcount + 1;
      ll = kk + length(ss)-1;
      pen.ki_val(kk:ll) = ss;
      pen.ki_row(kk:ll) = ii-1;
      pen.ki_col(kk:ll) = jj-1;
      kk = kk + length(ss);
   end
end


pen.ioptions = [1 50 100 2 0 1 0 0];
pen.foptions = [10000 0.7 0.1 1.0E-7 1.0E-6 1.0E-14 1.0E-3 1];