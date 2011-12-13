% function Prob = binbin2lin(Prob, idx4, idx1, idx2, idx3)
%
% Adds constraints when modeling with binary variables which is the product
% of two other variables.
%
% b4 = b1 * b2. The problem should be built with the extra variable b4 in
% place of the b1*b2 products. The indices of the unique product variables
% are needed to convert the problem properly.
%
% Three inequalities are added to the problem:
%
% b4 <= b1
% b4 <= b2
% b4 >= b1 + b2 - 1
%
% By adding this b4 will always be the product of b1 and b2.
%
% The routine also handles products of three binary variables.
%
% b4 = b1 * b2 * b3. The following constraints are then added:
%
% b4 <= b1
% b4 <= b2
% b4 <= b3
% b4 >= b1 + b2 + b3 - 1
%
% INPUT PARAMETERS
%
% Prob         Problem structure to be converted
% idx4         Indices for b4 variables
% idx1         Indices for b1 variables
% idx2         Indices for b2 variables
% idx3         Indices for b3 variables (optional)
%
% OUTPUT PARAMETERS
%
% Prob         Problem structure with added constraints

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2007 by Tomlab Optimization Inc., $Release: 5.7.0$
% Written Oct 4, 2005.   Last modified Mar 5, 2007.

function Prob = binbin2lin(Prob, idx4, idx1, idx2, idx3)

if nargin < 4
   error('binbin2lin requires at least 4 inputs');
end

if nargin < 5
   idx3 = [];
else
   if length(idx4) ~= length(idx3)
      error('idx3 and idx4 do not have the same length');
   end
   idx3 = full(idx3(:));
end

len4 = length(idx4);

if len4 ~= length(idx1)
   error('idx4 and idx1 do not have the same length');
end

if len4 ~= length(idx2)
   error('idx4 and idx2 do not have the same length');
end

idx1 = full(idx1(:));
idx2 = full(idx2(:));
idx4 = full(idx4(:));

A1 = sparse(repmat((1:len4)',2,1),[idx4;idx1],[ones(len4,1);-1*ones(len4,1)],len4,Prob.N);
A2 = sparse(repmat((1:len4)',2,1),[idx4;idx2],[ones(len4,1);-1*ones(len4,1)],len4,Prob.N);

if isempty(idx3)
   mAdd = len4*3;
   Prob.mLin = Prob.mLin + mAdd;
   A3 = sparse(repmat((1:len4)',3,1),[idx4;idx1;idx2],[-1*ones(len4,1);ones(len4,1);ones(len4,1)],len4,Prob.N);
   Prob.A = [Prob.A;A1;A2;A3];
   Prob.b_U = [Prob.b_U;zeros(len4*2,1);ones(len4,1)];
else
   mAdd = len4*4;
   Prob.mLin = Prob.mLin + mAdd;
   A3 = sparse(repmat((1:len4)',2,1),[idx4;idx3],[ones(len4,1);-1*ones(len4,1)],len4,Prob.N);
   A4 = sparse(repmat((1:len4)',4,1),[idx4;idx1;idx2;idx3],[-1*ones(len4,1);ones(len4,1);ones(len4,1);ones(len4,1)],len4,Prob.N);
   Prob.A = [Prob.A;A1;A2;A3;A4];
   Prob.b_U = [Prob.b_U;zeros(len4*3,1);2*ones(len4,1)];
end

Prob.b_L = [Prob.b_L;-inf*ones(mAdd,1)];

% MODIFICATION LOG
%
% 051004 med   Created
% 070305 med   Corrected constraints