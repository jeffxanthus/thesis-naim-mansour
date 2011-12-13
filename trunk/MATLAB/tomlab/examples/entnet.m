function [x,y,z,inform] = entnet(fname)
%        [x,y,z,inform] = entnet(fname)
% entnet loads network data from file 'fname'
% and solves a network problem with entropy objective:
%    min sum xj*ln(xj)  s.t.  Ax = b, x > 0.
%
% The network data in 'fname' is a list of  i j  pairs.
% The corresponding variables xij must satisfy Ax = b,
% where A and b are defined by
%    sum(i) xik  -  sum(j) xkj  = 0,
%                  sum(ij) xij  = 1.

% Michael Saunders, SOL, Stanford University.
%-----------------------------------------------------------------------
% 08 Feb 2002: First version of entnet.m
%              derived from John Tomlin's AMPL model.
%              Network data was built-in as an explicit matrix.
% 09 Feb 2002: entnet loads  i j  data from a file
%              and creates "A" as a sparse matrix.
%              Ideally, "A" should be an operator, but then pdsco.m
%              doesn't know how to do diagonal preconditioning for LSQR.
%              (We could make a special version of pdsco to do it.)
%-----------------------------------------------------------------------

A = load(fname);      % A(:,1) => +1,   A(:,2) => -1.
m = max(A(:,1));
k = min(A(:,1));
if k==0   % Adjust for 0 base.
   A = A + 1;   m = m+1;
end
n = size(A,1);        % Number of nodes
J = (1:n)';           % Column indices
m = m+1;              % Include extra constraint, sum(x) = 1
e = ones(n,1);
A = sparse([A(:,1); A(:,2); m*e], [J; J; J], [e; -e; e]);

% For entnet.dat, A will now be this:
% [ 1  1  0  0  0  0  0  0  0  0  0  0  0  0 -1
%  -1  0  1  1  0  0  0  0  0  0  0  0  0  0  0
%   0 -1  0  0  1  1  0  0  0  0  0  0  0  0  0
%   0  0 -1  0 -1  0  1  1  0  0  0  0  0  0  0
%   0  0  0 -1  0 -1  0  0  1  1  1  0  0  0  0
%   0  0  0  0  0  0 -1  0 -1  0  0  1  0  0  0
%   0  0  0  0  0  0  0 -1  0 -1  0  0  1  0  0
%   0  0  0  0  0  0  0  0  0  0 -1  0  0  1  0
%   0  0  0  0  0  0  0  0  0  0  0 -1 -1 -1  1
%   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1 ]

b     = zeros(m,1);   b(m)  = 1;
xsize = 1;            zsize = 1;
xsize = 0.2;          zsize = 0.01;     %%% TEST
x0    = e;            x0    = x0*xsize;
y0    = zeros(m,1);
z0    = e;            z0    = z0*zsize;

options        = pdscoSet;
options.x0min  = xsize;
options.z0min  = zsize;
options.mu0    = 1e-4;
%options.Method = 3; % Use SOL lsqr
options.Method = 4;  % Use Tomlab SOL lsqr MEX - Tlsqr
options.wait   = 1;
Prob.PriLevOpt = 1;
% To make pdsco quiet:
% options.wait   = 0;
% Prob.PriLevOpt = 0;
[x,y,z,inform] = pdsco('entropy',A,b,m,n,options,x0,y0,z0,xsize,zsize,Prob);
%keyboard

%-----------------------------------------------------------------------
% End function entnet
%-----------------------------------------------------------------------
