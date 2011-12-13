function snoptA_Ex2

%% snoptA_Ex2;

% Example of how to run a regular Tomlab problem in snoptA, but without exploiting
% linear variables.

Prob = probInit('con_prob',15);

G = ones(1,Prob.N);  % Put objective on first row
objRow = 1;

G = [ G ; estConsPattern(Prob,10) ]; % Nonlinear constraints

% Last rows are linear rows - there are no nonlinear elements there:
G = [ G ; spalloc(Prob.mLin,Prob.N,0) ];

% Linear coefficients for the 1+Prob.mNonLin first rows - NONE.
A = spalloc(1+Prob.mNonLin,Prob.N,0);
A = [ A ; sparse(Prob.A) ];

% Function bounds - objective, nonlinear, linear
b_L = [-Inf ; Prob.c_L(:) ; Prob.b_L(:) ];
b_U = [ Inf ; Prob.c_U(:) ; Prob.b_U(:) ];

x_L = Prob.x_L;
x_U = Prob.x_U;
x_0 = Prob.x_0;

Warmstart=0;
xState=[];
fState=[];
nS=[];

% optPar similar to snopt except: optPar(14:17) are not available in
% snoptA. Refer to snoptTL.m 
optPar = -999*ones(72,1);
% optPar(1) = 1;
% optPar(2) = 1;
optPar(13) = 1;

% Print and summary files, as snopt.
PrintFile='snpri.txt';
SummFile='snsum.txt';

% Print level: 0 - silent, 1 - summary, 2 - similar to print file 
PriLev=1;

% Callback file for solving regular TOMLAB problems with snoptA:
% Func = 'snoptATL_fg';
Func = @snoptATL_fg;


%% 
[x_k, xstate, xmul, F, Fmul, Fstate, Inform, nS, nInf, sInf] = ...
   snoptA( x_L, x_U, x_0, b_L, b_U, A, G, Func, objRow, optPar, ...
   Warmstart, xState, fState, nS, Prob, PrintFile, SummFile, PriLev);

          
          
