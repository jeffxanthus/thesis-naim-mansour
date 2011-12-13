%% example_sdp - tomSym SDP demonstration
% An SDP example, from Wikipedia article about semidefinite programming
%
% Minimize
%    rac
% subject to
%   [1 rab rac; rab 1 rbc; rac rbc 1] positive semidefinite
%   -0.2 <= rab <= -0.1
%    0.4 <= rbc <=  0.5

% Variables
toms rab rac rbc

% Objective function
f = rac;

% Constraints
M = [1 rab rac; rab 1 rbc; rac rbc 1]
con = { positiveSemidefinite(M), ...
    -0.2 <= rab <= -0.1, ...
    0.4 <= rbc <=  0.5};

% Initial conditions
x0 = struct('rab',-0.15,'rbc',0.45,'rac',0);

%% Compile and solve problem as a NLP
options = struct;
options.type = 'con';
options.name = 'SDP example';
solution = ezsolve(f,con,x0,options);

% Evaluate M using the returned solution
M_opt = subs(M,solution)

%% Compile and solve problem as a SDP
options.type = 'sdp';
solution = ezsolve(f,con,x0,options);

% Evaluate M using the returned solution
M_opt = subs(M,solution)

%% Check the eigenvalues
% The eigenvalues should all be positive, even though small negative ones
% may occur due to tolerances.
eigenvalues = eigs(M_opt)