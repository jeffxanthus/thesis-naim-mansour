% jogo09_prob.m
%
% Unconstrained and constrained global optimization problems from 
% "Handling Box, Linear and Quadratic-Convex Constraints for 
% Boundary Optimization with DE and DIRECT Algorithms"
% by Massimo Spadoni and Luciano Stefanini, University of Urbino, Italy
%
%
% function [probList, Prob] = jogo09_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%    n      Defines the number of variables for problems 4 to 8 and 14 to 25.
%           Must be an even number for problem 5.
%           Default values are used if n is not given as input.
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 9, 2009.   Last modified July 8, 2009.

function [probList, Prob] = jogo09_prob(P, varargin)

if nargin < 1
    P=[];
end


probList=str2mat(...
     'Low-dimensional 1 - Mladenovic'...
    ,'Low-dimensional 2 - Sun'...
    ,'Low-dimensional 3 - Churilov'...
    ,'Box constrained with boundary solutions 1 - Cox'...
    ,'Box constrained with boundary solutions 2 - Maranas - Circle Packing'...
    ,'Box constrained with boundary solutions 3 - Locatelli 1'...
    ,'Box constrained with boundary solutions 3 - Locatelli 2'...
    ,'Box constrained with boundary solutions 3 - Locatelli 3'...
    ,'Fully constrained with boundary solutions 1'...
    ,'Fully constrained with boundary solutions 2'...
    ,'Fully constrained with boundary solutions 3'...
    ,'Fully constrained with boundary solutions 4'...
    ,'Fully constrained with boundary solutions 5'...
    ,'Fully constrained with boundary solutions 6'...
    ,'Fully constrained with boundary solutions 7'...
    ,'Fully constrained with boundary solutions 8 - Random constraints'...
    ,'Fully constrained with boundary solutions 9 - Random constraints'...
    ,'Fully constrained with boundary solutions 10 - Random constraints'...
    ,'Fully constrained with boundary solutions 11 - Random constraints'...
    ,'Fully constrained with boundary solutions 12 - Random constraints'...
    ,'Fully constrained with boundary solutions 13 - Random constraints'...
    ,'Fully constrained with boundary solutions 14 - Random constraints'...
    ,'Fully constrained with boundary solutions 15 - Random constraints'...
    ,'Fully constrained with boundary solutions 16 - Random constraints and objective'...
    ,'Fully constrained with boundary solutions 17 - Random constraints and objective'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return;
end

n = varargin{2};
if isempty(n) | n <= 0
  % number of variables to use in problem 4
  switch P
     case {4,6,7,8,14,15}
        n = 5;
     case {4}
        n = 5;
     case {5}
        n = 10; % 2*5
     otherwise
        n = 10;
  end
end
n

IntVars = []; f_Low = [];
x_0  = []; % Currently not used
user = [];
if P == 1
    Name = 'Low-dimensional 1 - Mladenovic';
    x_L = [0  0 ]';
    x_U = [6 Inf]';
    b_L = [ 0  -6 ]';
    b_U = [];
    A = [1/sqrt(3) -1;
         -1  -sqrt(3)];
    c_L = []; c_U = [];
    x_opt = [];
    f_opt = -1;
    x_min = x_L;
    x_max = x_U;
elseif P == 2
    Name = 'Low-dimensional 2 - Sun';
    x_L = [2 2]';
    x_U = [6 6]';
    b_L = [-44 -48]'; 
    b_U = []; 
    A = [-2 -8;
         -7 -4];
    c_L = []; c_U = [];
    x_opt = [];
    f_opt = -200.1667;
    x_min = [2 2];
    x_max = [6 6];
elseif P == 3
    Name = 'Low-dimensional 3 - Churilov';
    x_L = [0 0 0]';
    x_U = [Inf Inf Inf]';
    b_L = [-1 -9 -12 -3 -3.5 7]'; b_U = [];
    A = [ 2 -1   -1;
          0 -3   -1;
         -3 -1.5 -1;
         -3  3   -1;
         -1  1   -2;
          2  3   -1];
    c_L = []; c_U = [];
    x_opt = [];
    f_opt = 6.5086;
    x_min = x_L;
    x_max = x_U;
elseif P == 4
    Name = 'Box constrained with boundary solutions 1 - Cox - Random noise';
    x_L = -2*ones(n,1);
    x_U = 2*ones(n,1);
    b_L = []; b_U = []; A = [];
    c_L = []; c_U = [];
    x_opt = [];
    switch n
      case 5
         f_opt = -87.3898;
      case 10
         f_opt = -169.6234;
      otherwise
         f_opt = [];
    end
    x_min = x_L;
    x_max = x_U;
elseif P == 5
    Name = 'Box constrained with boundary solutions 2 - Maranas - Packing circles in a square';
    if mod(n,2) ~= 0
       error('n must be an even number (variables correspond to the center coordinates of circles)');
    end
    x_L = zeros(n,1);
    x_U = ones(n,1);
    b_L = []; b_U = []; A = [];
    c_L = []; c_U = [];
    x_opt = [];
    switch n
      case 5
         f_opt = -0.7071;
      case 10
         f_opt = -0.5;
      otherwise
         f_opt = [];
    end
    x_min = x_L;
    x_max = x_U;
elseif P == 6
    Name = 'Box constrained with boundary solutions 3 - Locatelli 1';
    x_L = zeros(n,1);
    x_U = ones(n,1);
    b_L = []; b_U = []; A = [];
    c_L = []; c_U = [];
    x_opt = [];
    switch n
      case 5
         f_opt = -5.1962;
      case 10
         f_opt = -6.7082;
      otherwise
         f_opt = [];
    end
    x_min = x_L;
    x_max = x_U;
elseif P == 7
    Name = 'Box constrained with boundary solutions 3 - Locatelli 2';
    x_L = zeros(n,1);
    x_U = ones(n,1);
    b_L = []; b_U = []; A = [];
    c_L = []; c_U = [];
    x_opt = [];
    switch n
      case 5
         f_opt = -2.7146;
      case 10
         f_opt = -4.0079;
      otherwise
         f_opt = [];
    end
    x_min = x_L;
    x_max = x_U;
elseif P == 8
    Name = 'Box constrained with boundary solutions 3 - Locatelli 3';
    x_L = zeros(n,1);
    x_U = ones(n,1);
    b_L = []; b_U = []; A = [];
    c_L = []; c_U = [];
    x_opt = [];
    switch n
      case 5
         f_opt = -2.6259;
      case 10
         f_opt = -4.5096;
      otherwise
         f_opt = [];
    end
    x_min = x_L;
    x_max = x_U;
elseif P == 9
    Name = 'Fully constrained with boundary solutions 1';
    x_L = [0 0 0 0 0 0 0 0 0   0   0   0 0]';
    x_U = [1 1 1 1 1 1 1 1 1 100 100 100 1]';
    b_L = [-10 -10 -10 0 0 0 0 0 0]';
    b_U = [];
    A = [-2 -2  0 0 0 0 0 0 0 -1 -1  0 0;
         -2  0 -2 0 0 0 0 0 0 -1  0 -1 0;
          0 -2 -2 0 0 0 0 0 0  0 -1 -1 0;
          0  0  0 2 1 0 0 0 0 -1  0  0 0;
          0  0  0 0 0 2 1 0 0  0 -1  0 0;
          0  0  0 0 0 0 0 2 1  0  0 -1 0;
          8  0  0 0 0 0 0 0 0 -1  0  0 0;
          0  8  0 0 0 0 0 0 0  0 -1  0 0;
          0  0  8 0 0 0 0 0 0  0  0 -1 0];
    c_L = []; c_U = [];
    x_opt = [1 1 1 1 1 1 1 1 1 3 3 3 1];
    f_opt = -15;
    x_min = x_L;
    x_max = x_U;
elseif P == 10
    Name = 'Fully constrained with boundary solutions 2';
    user.Q = [1 0 1 1 0;
              0 1 0 1 1;
              1 0 1 0 1;
              1 1 0 1 0;
              0 1 1 0 1];
    x_L = [0 0 0 0 0]';
    x_U = [1 1 1 1 1]';
    b_L = [1]; 
    b_U = []; 
    A = [1 1 1 1 1];
    c_L = []; c_U = [];
    x_opt = [];
    f_opt = 0.5;
    x_min = x_L;
    x_max = x_U;
elseif P == 11
    Name = 'Fully constrained with boundary solutions 3';
    user.Q = [1 0 0 0 0 0 1 1 1 1 1 1
              0 1 0 0 1 1 0 0 1 1 1 1
              0 0 1 1 0 1 0 1 0 1 1 1
              0 0 1 1 1 0 1 0 1 0 1 1
              0 1 0 1 1 0 1 1 0 1 0 1
              0 1 1 0 0 1 1 1 1 0 0 1
              1 0 0 1 1 1 1 0 0 1 1 0
              1 0 1 0 1 1 0 1 1 0 1 0
              1 1 0 1 0 1 0 1 1 1 0 0
              1 1 1 0 1 0 1 0 1 1 0 0
              1 1 1 1 0 0 1 1 0 0 1 0
              1 1 1 1 1 1 0 0 0 0 0 1];
    x_L = [0 0 0 0 0 0 0 0 0 0 0 0]';
    x_U = [1 1 1 1 1 1 1 1 1 1 1 1]';
    b_L = [1]; 
    b_U = []; 
    A = [1 1 1 1 1 1 1 1 1 1 1 1];
    c_L = []; c_U = [];
    x_opt = [];
    f_opt = 0.333333333;
    x_min = x_L;
    x_max = x_U;
elseif P == 12
    Name = 'Fully constrained with boundary solutions 4';
    user.Q = [14 15   16    0    0
              15 14   12.5 22.5 15
              16 12.5 10   26.5 16
               0 22.5 26.5  0    0
               0 15   16    0   14];
    x_L = [0 0 0 0 0]';
    x_U = [1 1 1 1 1]';
    b_L = [1]; 
    b_U = [];
    A = [1 1 1 1 1];
    c_L = []; c_U = [];
    x_opt = [];
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 13
    Name = 'Fully constrained with boundary solutions 5';
    x_L = [-10 -10 -10 -10 -10 -10 -10 -10 -10 -10]';
    x_U = [10 10 10 10 10 10 10 10 10 10]';
    b_L = []; 
    b_U = [105 0 12]'; 
    A = [ 4  5 0 0 0 0  -3 9 0  0;
         10 -8 0 0 0 0 -17 2 0  0;
         -8  2 0 0 0 0   0 0 5 -2];
    c_L = [];
    c_U = [0 0 0 0 0];
    % note that x_opt somewhat violates b_U(1), b_U(3) and c_U(1:3)
    x_opt = [2.172 2.364 8.774 5.096 0.991 1.431 1.322 9.829 8.280 8.376];
    f_opt = 24.306;
    x_min = x_L;
    x_max = x_U;
elseif P == 14
    Name = 'Fully constrained with boundary solutions 6';
    x_L = zeros(n,1);
    x_U = ones(n,1);
    b_L = []; b_U = []; A = [];
    c_L = [];
    c_U = [1];
    x_opt = 1/sqrt(n)*ones(n,1);
    switch n
      case 5
         f_opt = -1.0; %1.0 in paper, they seem to have excluded the minus sign.
      case 10
         f_opt = -1.0; %1.0 in paper, they seem to have excluded the minus sign.
      otherwise
         f_opt = [];
    end
    x_min = x_L;
    x_max = x_U;
elseif P == 15
    Name = 'Fully constrained with boundary solutions 7';
    x_L = zeros(n,1);
    x_U = ones(n,1);
    b_L = []; b_U = [1];
    A = ones(1,n);
    c_L = [];
    c_U = [];
    x_opt = 1/n*ones(1,n);
    switch n
      case 5
         f_opt = -1.0; %1.0 in paper, they seem to have excluded the minus sign.
      case 10
         f_opt = -1.0; %1.0 in paper, they seem to have excluded the minus sign.
      otherwise
         f_opt = [];
    end
    x_min = x_L;
    x_max = x_U;
elseif any(P == [16 17 18 19 20 21 22 23 24 25])
    Name = ['Fully constrained with boundary solutions ' num2str(P-8) ' - Random constraints'];
    switch P % Define random parameters for objective function
       case 24
          Name = [Name ' and objective'];
          user.beta = rand([n,1]);
          user.j1 = ceil(n.*rand([n,1]));
          user.j2 = ceil(n.*rand([n,1]));
       case 25
          Name = [Name ' and objective'];
          k = 2; % Must be >= 1. k = 2 is used in paper
          z = sort([0;k*rand([n-1,1]);k]);
          user.alfa = z(2:n+1)-z(1:n);
    end
    x_L = zeros(n,1);
    x_U = 100*ones(n,1);
if n~=10
    mL = 25; % number of randomly generated linear constraints
    % bounds on random numbers
    alfa1 = 1; alfa2 = 11;
    beta1 = -10; beta2 = 10;
    gamma1 = 2; gamma2 = 22;
    % In the paper:
    % q_i = -sum(p_ij) + [rand(gamma1,gamma2)]     for j = 1...mL
    %   where -p_1j belongs to [rand(alfa1,alfa2)]
    %     and -p_ij belongs to [rand(beta1,beta2)] for i = 2...n
    %
    % This is inconsistent with the values for A and b_L presented on page 18.
    % Instead we use 
    % q_i = sum(p_ij) - [rand(gamma1,gamma2)]
    %
    % q_i is the lower bound on the constraint i
    A = -floor([ alfa1+(alfa2-alfa1).*rand([1,   n]) ;
                 beta1+(beta2-beta1).*rand([mL-1,n])  ]); % generated randomly
    b_L = [ sum(A,2)-floor(gamma1+(gamma2-gamma1).*rand([mL,1]))]; % generated randomly
    b_U = [];
    f_opt = [];
else % n=10, predetermined matrix
    %x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 RHS
    Z=[...
    -10 -3 -7 -5 -9 -8 -5 -1 -9 -5 -70
    -2 1 -7 0 -1 6 10 1 4 4 14
    -5 2 -1 -7 -5 -2 3 -1 3 7 -16
    -8 -6 1 7 0 4 -3 8 3 7 -3
    -4 0 -7 -9 -2 -9 9 1 -1 7 -31
    7 6 -6 5 6 -4 10 -4 8 2 10
    2 -3 -2 5 3 2 -2 -7 10 -7 -17
    -8 -6 -6 -7 -5 -4 -2 5 1 1 -39
    -8 10 -3 -4 -3 5 10 5 -7 -6 -15
    2 -3 4 8 1 2 10 -7 -8 1 -6
    -7 3 5 10 -1 -8 7 6 5 1 17
    9 -6 4 -7 -5 -3 -1 -6 7 1 -27
    3 0 0 7 9 6 9 -8 -7 2 9
    -6 -4 -4 5 -2 -6 3 6 6 -8 -24
    10 2 4 -3 9 -2 -2 6 -2 10 26
    8 4 -6 5 2 8 -4 10 -9 5 11
    6 7 -1 1 4 6 -3 9 -3 10 16
    7 7 3 9 -7 -2 9 -2 -7 -3 6
    -2 -3 -4 -9 10 -2 1 7 10 -3 -9
    5 4 0 -1 -5 3 2 -6 8 -9 -7
    7 0 2 2 -9 -1 3 7 -6 -1 -10
    10 7 -3 0 -9 1 7 7 2 2 10
    -4 -3 -2 4 -5 10 -3 -9 -7 7 -20
    2 3 -5 2 2 10 -3 2 -4 -2 -3
    -8 -7 -9 6 1 4 -4 4 -3 -4 -34
    ];
    A=Z(:,1:10);
    b_U = [];
    b_L = Z(:,11);
    switch P
      case 16
        f_opt = -41.7017;
      case 17
        f_opt =  -9.7248;
      case 18
        f_opt = -107.2079;
      case 19
        f_opt = -169.3058;
      case 20
        f_opt = -4.9099;
      case 21
        f_opt = 0.3255;
      case 22
        f_opt = -203.4805;
      case 23
        f_opt = -41.1397;
      case 24
        f_opt = -894.8621;
      case 25
        f_opt = -2.1334;
    end
end
    c_L = [];
    c_U = [];
    x_opt = [];
    x_min = x_L;
    x_max = x_U;
    x_0 = ones(n,1); % should always be feasible
else
    error('jogo09_prob: Illegal problem number');
end

% Set x_0 to zeros
if isempty(x_0)
    x_0=zeros(length(x_L),1);
end

% Define the Prob
c  = 'jogo09_c';
if isempty(c_L) & isempty(c_U)
    c  = [];
end

Prob = glcAssign('jogo09_f', x_L, x_U, Name, A, b_L, b_U, ...
    c, c_L, c_U, x_0, IntVars, [], [], [], ...
    f_Low, x_min, x_max, f_opt, x_opt);
Prob.P = P;
Prob.user = user;

% MODIFICATION LOG
%
% 090609  bjo  First version created
% 090708  hkh  Added f_opt values, fixed A,b_L for n=10, P=16:25
