% function Prob = bincont2lin(Prob, idx_prod, idx_bin, idx_cont)
%
% Adds constraints when modeling with binary variables which are multiplied
% by integer or continuous variables. This is the most efficient way to
% get rid off quadratic objectives or constraints.
%
% prod = bin * cont. The problem should be built with the extra variables
% prod in place of the bin*cont products. The indices of the unique product
% variables are needed to convert the problem properly.
%
% Three inequalities are added to the problem:
%
% prod <= cont
% prod >= cont - x_U * (1 - bin)
% prod <= x_U * bin
% prod >= x_L * bin
%
% By adding this prod will always equal bin * cont.
%
% INPUT PARAMETERS
%
% Prob         Problem structure to be converted
% idx_prod     Indices for product variables
% idx_bin      Indices for binary variables
% idx_cont     Indices for continuous/integer variables
%
% OUTPUT PARAMETERS
%
% Prob         Problem structure with added constraints

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2008 by Tomlab Optimization Inc., $Release: 6.1.0$
% Written Oct 4, 2005.   Last modified May 22, 2008.

function Prob = bincont2lin(Prob, idx_prod, idx_bin, idx_cont)

if nargin < 4
    error('bincont2lin requires 4 inputs');
end

len = length(idx_prod);

if len ~= length(idx_bin)
    error('idx_prod and idx_bin do not have the same length');
end

if len ~= length(idx_cont)
    error('idx_prod and idx_cont do not have the same length');
end

if ~any(isfinite(Prob.x_U(idx_cont)))
    error('An upper bound is required on the integer/continuous variables');
end

if any(abs(Prob.x_U(idx_cont) - Prob.x_U(idx_prod)) > 1e-6)
    error('The upper bounds need to be same for the integer/continuous and products');
end

idx_prod = full(idx_prod(:));
idx_bin  = full(idx_bin(:));
idx_cont = full(idx_cont(:));

mAdd = len*4;
Prob.mLin = Prob.mLin + mAdd;

x_U = Prob.x_U(idx_cont);
x_L = Prob.x_L(idx_cont);

A1 = sparse(repmat((1:len)',2,1),[idx_prod;idx_cont],[ones(len,1);-1*ones(len,1)],len,Prob.N);
A2 = sparse(repmat((1:len)',3,1),[idx_cont;idx_bin;idx_prod],[ones(len,1);x_U;-1*ones(len,1)],len,Prob.N);
A3 = sparse(repmat((1:len)',2,1),[idx_prod;idx_bin],[ones(len,1);-x_U],len,Prob.N);
A4 = sparse(repmat((1:len)',2,1),[idx_prod;idx_bin],[ones(len,1);-x_L],len,Prob.N);

bU = zeros(mAdd,1);
bU(len+1:len*2) = x_U;
bU(3*len+1:end) = inf*ones(len,1);
bL = [-inf*ones(mAdd-len,1);zeros(len,1)];
Prob.A = [Prob.A;A1;A2;A3;A4];
Prob.b_U = [Prob.b_U;bU];
Prob.b_L = [Prob.b_L;bL];

% MODIFICATION LOG
%
% 051004 med  Created.
% 070220 med  Updated b_U vector
% 080522 med  Lower bounds now supported as well 