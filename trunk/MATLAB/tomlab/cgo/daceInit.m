% daceInit.m
%
% daceInit finds an initial grid using a Latin Hypercube space-filling design.
%
% function X = daceInit(M, s, x_L, x_U, X, minDist);
%
% INPUT:
%
% M       Number of points if M > 0. 
%         If M == [] M = d+1, where d = dim(x_L)
%         If M == 0  M = (d+1)*(d+2)/2, where d = dim(x_L)
%         If M <  0, k = abs(M)
%         k is in principle the dimension. It determines how many points to use
%         k    # of points
%         1    21
%         2    21
%         3    33
%         4    41
%         5    51
%         6    65
%        >6    65
%
% s       How many intervals s to create in each dimension
%         Default s = min(100,M)
%
% x_L     Lower bounds for each element in x.
% x_U     Upper bounds for each element in x.
%
% X       If X is given, the M new points are added to X,
%         except in the case:  k==2 & d==2 & Method == 1
%
% MinDist If not [], all new points must have |xNew - X(:,i)| >= minDist
%
% M         Number of points to generate, if to overrule the k value
%           Default empty, i.e. M is not used, only k 
%
% OUTPUT:
% X         Sample point matrix, d x M

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Oct 29, 1998.   Last modified Sept 17, 2009.

function X = daceInit(M, s, x_L, x_U, X, minDist, Method)

if nargin < 7 
   Method = []; 
   if nargin < 6 
      minDist = [];
      if nargin < 5 
         X = [];
      end
   end
end

if isempty(Method), Method = 2; end

d = length(x_L);
k = 0;

if isempty(M) 
   M = d+1; 
elseif M < 0
   k = abs(M);
   switch k
     case {1,2}
      M = 21;
     case 3
      M = 33;
     case 4
      M = 41;
     case 5
      M = 51;
     case 6
      M = 65;
     otherwise
      M = 65;
   end
elseif M == 0
   M = (d+1)*(d+2)/2;
end
if isempty(s), s = min(100,M); end

if k==2 & d==2 & Method == 1
   n = 21;           % number of sampled points
   %u = 0:0.05:1;
   v = [0.5 0.15 0.75 1 0.35 0.6 0.1 0.85 0.3 0.55 0.05 ...
        0.95 0.7 0.25 0.45 0 0.9 0.65 0.2 0.4 0.8];
   X = zeros(2,n);
   X(1,:) = x_L(1) + ( x_U(1) - x_L(1) )*(0:0.05:1);
   X(2,:) = x_L(2) + ( x_U(2) - x_L(2) )*v;
   return
end

xD = x_U-x_L;
   
if isempty(X)
   X  = zeros(d,M);
   k  = 0;
else
   k  = size(X,2);
   X  = [X,zeros(d,M)];
end

%rand('state',2);
if minDist > 0
   % Generate M Latin Hypercube points based on sampling blocks of s 
   % random points in d-dimensional space
   % Each block of s points from one of s subdivisions of [0,1], 
   % all with uniform distribution 
   % All M points must be at least the distance minDist from the given set X
   ok  = 0;
   pnt = 0;
   T   = zeros(d,s);
   while ok < M & pnt < 100000
       for i = 1:d
           T(i,:) = x_L(i) + xD(i)*((rand(1,s)+(randperm(s)-1))/s);
       end
       for i = 1:s
           xNew = T(:,i);
           %NHQ Now possible to use minDist with an empty X
           if k == 0
              % No point sampled yet, accept point
              ok = ok + 1;
              X(:,k+ok) = xNew;
           elseif min(tomsol(30,xNew,X(:,1:k))) >= minDist
              % Distance to every X > minDist, accept point
              ok = ok + 1;
              X(:,k+ok) = xNew;
           end
           if ok >= M, break; end
       end
       pnt = pnt + s;
   end
   if pnt >= 100000
      X = X(:,1:k+ok);
   end
elseif Method == 1
   % Generate M Latin Hypercube points based on s permutated fixed steps 
   % in d-dimensional space
   l = k;
   %n = ceil(M/s);
   v = linspace(0,1,s);
   for j = 1:ceil(M/s)
       for i = 1:d
           [a b]=sort(rand(1,s));   
           X(i,l+1:l+s) = x_L(i) + xD(i)*v(b);
       end
       l = l + s;
   end
elseif Method == 2
   % Generate M Latin Hypercube points based on sampling blocks of s 
   % random points in d-dimensional space
   % Each block of s points from one of s subdivisions of [0,1], 
   % all with uniform distribution 
   l = k;
   for j = 1:ceil(M/s)
       for i = 1:d
           X(i,l+1:l+s) = x_L(i) + xD(i)*((rand(1,s)+(randperm(s)-1))/s);
       end
       l = l + s;
   end
end

%if k == 2  % Plot initial points
%   plot(X(1,:),X(2,:),'*r');
%   pause
%end

% MODIFICATION LOG:
%
% 980626  hkh  Avoid feval
% 010715  hkh  Only compute X, use lower and upper bounds as input
% 020831  hkh  Use rand('state'), Matlab 5.x random generator
% 040306  hkh  Avoid setting the state, done in ego and rbfSolve 
% 040307  hkh  Add extra parameter M, to generate a large set of points
% 050117  med  mlint revision
% 050322  hkh  Wrong index accessed for n=2, x_L,x_U index 1 instead of 2
% 050427  hkh  Complete revision, correction of algorithm
% 050427  hkh  New Method = 2 is better to use
% 050427  hkh  Generate M points with distance > minDist to set X
% 050509  hkh  Sign error in new method, generated negative numbers
% 080410  hkh  Revised. Now possible to use minDist with an empty X
% 090821  hkh  Remove unnecessary init of y = zeros(n,1)
% 090901  hkh  Correct comments, input X not used in special case
% 090917  hkh  Limit default s, if isempty(s), s = min(100,M);
