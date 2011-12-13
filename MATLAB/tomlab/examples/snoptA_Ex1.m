function snoptA_Ex1

% snoptA_Ex1.m
%
% min f1 = exp(x(2))*x(4) + x(2)^2 + 3*x(1)-x(3)+4*x(4)+x(5)
%  x 
% 
% s/t    f2 == 2 (equality)
%        f3 >= 0
% 
% where f2 = x(3)^2 + sin(x(4)) + x(2) - 3*x(5) == 2
%       f3 = 0 + x(1)-x(3)
%

% Linear coefficient matrix. There are two variables that occur only in
% linear expressions in f1,f2,f3; those are x(1) and x(5).
A = [ 3  0 -1  4  1  ; ...
      0  1  0  0 -3  ; ...
      1  0 -1  0  0 ];

A = sparse(A);
[nF,n] = size(A);

% G is similar to ConsPattern in normal TOMLAB NLP problems, but 
% also includes a row for the objective function. 
i = [1,2,1,2];
j = [2,3,4,4];
v = 1;
G = sparse(i,j,v,nF,n);

% NOTE: A and G should be of identical size (unless one or the other is
% empty). They overlap to form the Jacobian for the entire problem:
%
% F(x) = f(x) + A*x
% J(x) = G(x) + A
%
% so that G is the nonlinear terms of the Jacobian (or the Jacobian of
% f(x)) and A is the linear terms in all functions. It is important to code
% your problem so that you do not accidentally repeat a linear term in both
% G and A. 

% Variable bounds. Should be nonempty, equal-sized, and match the column
% dimension of A and G.
xl = -5*ones(n,1);
xu =  5*ones(n,1);

% Function bound vectors. 
%
% Unless the problem is one of finding a feasible point,
% one of the rows represent the objective function, indicated by objRow.
% This row should be unbounded in bl,bu. Dimension should be equal to the
% number of rows in A and G. 

bl = [-inf, 2, 0   ]';
bu = [ inf, 2, inf ]';
objRow = 1;

% optPar similar to snopt except: optPar(14:17) are not available in
% snoptA. Refer to snoptTL.m 
optPar = -999*ones(72,1);
optPar(1) = 1;
optPar(2) = 1;
optPar(13) = 3;

% Warmstart inputs (see below for example of warmstart)
Warm   = 0;
xstate = [];
fstate = [];
nS     = [];

% Prob - not used by snoptA mex, can be anything right now. 
Prob = struct('P',1);

% Starting point can be given
x0=[];

% Print and summary files, as snopt.
PrintFile='snprint.txt';
SummFile='snsum.txt';

% Print level: 0 - silent, 1 - summary, 2 - similar to print file 
PriLev=2;

% Function string (handle also allowed)
Func = 'snoptA_Ex1_fg';


[x, xstate, xmul, F, Fmul, Fstate, Inform, nS, nInf, sInf] = ...
   snoptA( xl, xu, x0, bl, bu, A, G, Func, objRow, optPar, ...
   Warm, xstate, fstate, nS, Prob, PrintFile, SummFile, PriLev);

%% Warmstart - uncomment to try this:
% try
%    x0 = x;
%    Warm = 1;
%    PrintFile='snpri.warm.txt';
%    SummFile='snsum.warm.txt';
%    [x,xs,xm,F,Fm,Fs,Inform,nS,nInf,sInf]=...
%       snoptA(xl,xu,x0,bl,bu,A,G,Func,objRow,optPar, ...
%       Warm,xstate,Fstate,nS,Prob,PrintFile,SummFile);
% catch
%    % speeds up development process to have everything cleared for rebuild...
%    lasterr
%    clear mex
% end
