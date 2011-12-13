%%*****************************************************************************
%% PENFEAS checks feasibility of a system of LMIs
%%
%% [ifeas,feas,xfeas] = penfeas_sdp(pen, options)
%%
%% Input: pen...system of LMIs in PENSDP format
%%        options(1)...output level of PENSDP (default = 0)
%%        options(2)...upper bound on the infty-norm of x (default = 1000)
%%                     needed when the system is unbounded
%%                     when set < 0 , no upper bounds are applied
%%
%% Output: ifeas...feasibility of the system: 0...system strictly feasible
%%                                            1...system (not strictly) feasible
%%                                           -1...system probably infeasible
%%
%%         feas...value of the minimized maximal eigenvalue of LMI
%%         xfeas...feasible point
%%
%% Copyright (c) 2003 by M. Kocvara and M. Stingl, PENOPT GbR
%% Version 10/01/2003
%%
%% Adapted to Tomlab /PENSDP: Anders Goran, Tomlab Optimization Inc.
%% Last modified Nov 7, 2003
%%
%%*****************************************************************************
function [ifeas,feas,xfeas] = penfeas_sdp(PEN,options)

if nargin < 1, error('Input arguments missing'); end  

if nargin < 2, 
   options(1) = 0;        %output level
   options(2) = 1000;     %box counstraints
end

p=penfeas1(PEN,options(2));

outlev= options(1);
p.ioptions = [1 50 100 outlev 0 0 0 0];
p.foptions = [100 0.7 0.1 1.0E-7 1.0E-6 1.0E-14 1.0E-2 1.0];

ifeas = -1;
[xf,feas] = pen(p);
if feas < 0, ifeas = 0; end
if abs(feas) < 1.0e-6 , ifeas = 1; end

xfeas = xf(1:length(xf)-1);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pen=penfeas1(PEN,upp_b)

pen.vars    = PEN.vars+1;
pen.constr  = PEN.constr;
pen.mconstr = PEN.mconstr;
pen.msizes  = PEN.msizes;

pen.x0  = zeros(1,PEN.vars+1);
pen.x0(PEN.vars+1)=-1;

pen.fobj = [zeros(1,PEN.vars),1];

%pen.ci= PEN.ci;
%pen.bi_dim = PEN.bi_dim;
%pen.bi_idx = PEN.bi_idx;
%pen.bi_val = PEN.bi_val;

%forget the original linear constraints, add box constraints
if upp_b < 0 
   pen.constr = 0;
   pen.bi_dim = 0;
   pen.bi_ci = 0;
   pen.bi_idx = 0;
   pen.bi_val = 0;
else   
   pen.constr = 2*PEN.vars;
   pen.ci= [0];
   pen.ci= upp_b*ones(1,2*PEN.vars);  % VALUE OF +- BOUNDS ON VARIABLES
   pen.bi_dim=ones(1,2*PEN.vars);
   for i=1:PEN.vars
      pen.bi_idx(1,2*i-1) = i-1;
      pen.bi_val(1,2*i-1) =-1;
      pen.bi_idx(1,2*i)   = i-1;
      pen.bi_val(1,2*i)   = 1; 
   end
end
% end of "forget"

PEN.ai_idx=reshape(PEN.ai_idx,1,length(PEN.ai_idx));
PEN.ai_nzs=reshape(PEN.ai_nzs,1,length(PEN.ai_nzs));
PEN.ai_val=reshape(PEN.ai_val,1,length(PEN.ai_val));
PEN.ai_col=reshape(PEN.ai_col,1,length(PEN.ai_col));
PEN.ai_row=reshape(PEN.ai_row,1,length(PEN.ai_row));

pen.ai_dim=PEN.ai_dim+1;
k=1;kk=1;
for i=1:PEN.mconstr
   pen.ai_idx(kk:kk+PEN.ai_dim(i))=[PEN.ai_idx(k:k+PEN.ai_dim(i)-1),PEN.vars+1];
   k = k + PEN.ai_dim(i);
   kk = k + i;
end   

gs = sum(PEN.msizes);

k=1;kk=1;
for i=1:PEN.mconstr
   pen.ai_nzs(kk:kk+PEN.ai_dim(i))=[PEN.ai_nzs(k:k+PEN.ai_dim(i)-1),PEN.msizes(i)];
   k = k + PEN.ai_dim(i);
   kk = k + i;
end

pen.ai_val = zeros(2,1);
pen.ai_col = zeros(2,1);
pen.ai_row = zeros(2,1);

k=1;kk=1;ksu=1;kksu=1;
for i=1:PEN.mconstr
   LL = sum(PEN.ai_nzs(k:k+PEN.ai_dim(i)-1));
   LL1 = sum(pen.ai_nzs(kk:kk+pen.ai_dim(i)-1));
   pen.ai_val(kksu:kksu+LL1-1)=[PEN.ai_val(ksu:ksu+LL-1)';-100*ones(PEN.msizes(i),1)];
   pen.ai_col(kksu:kksu+LL1-1)=[PEN.ai_col(ksu:ksu+LL-1)';(0:PEN.msizes(i)-1)'];
   pen.ai_row(kksu:kksu+LL1-1)=[PEN.ai_row(ksu:ksu+LL-1)';(0:PEN.msizes(i)-1)'];
   k = k + PEN.ai_dim(i);
   kk = k + i;
   ksu = ksu + LL;
   kksu = kksu + LL1;
end
