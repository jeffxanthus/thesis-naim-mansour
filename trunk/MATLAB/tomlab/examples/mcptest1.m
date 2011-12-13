% mcptest1
% Transportation model as a variational inequality in MATLAB with isoelstic
% demand.
% 
% lcptest1 is similar, see problem formualtion in this file.
%
% 0 <= x_k _|_ F(x) >= 0  (_|_ means complements)
%
% s.t b_L <= A <= b_U

clear all

MPEC = [1 0 0 0 1 0; ...
        2 0 0 0 2 0; 
        3 0 0 0 3 0; 
        4 0 0 0 4 0; 
        5 0 0 0 5 0; 
        6 0 0 0 6 0; 
        7 0 0 0 7 0; 
        8 0 0 0 8 0; 
        9 0 0 0 9 0; 
        10 0 0 0 10 0;         
        11 0 0 0 11 0];

d2c = [];
ConsPattern = [];

c_L = zeros(11,1);
c_U = inf*ones(11,1);

Prob = mcpAssign([], [], [], [], [], [], 'MCP 1', ones(11,1), ...
                          MPEC, [], ...
                          [], [], [], 'mcptest1_c', 'mcptest1_dc', d2c, ConsPattern, c_L, c_U);

Prob.KNITRO.options.ALG = 1;
R = tomRun('knitro', Prob, 1);                      
                      
% function Prob = mcpAssign(F, J, JacPattern, x_L, x_U, Name, x_0, ...
%                           A, b_L, b_U, x_min, x_max, f_opt, x_opt);

% Prob = mcpAssign('mcptest1_f', [], [], [], [], 'MCP 1', ones(11,1));
% Prob.PriLevOpt = 1;
% 
% Result = tomRun('path',Prob,1);
% pause

% -------------------------------------------------------------------
% Also give the Jacobian explicitly.
% 
% Prob2 = mcpAssign('mcptest1_f', 'mcptest1_J', [], [], [], 'MCP 2', ones(11,1));
% Prob2.PriLevOpt = 1;
% 
% Result2 = tomRun('path',Prob2,1);
% pause

% -------------------------------------------------------------------
% Test first problem with MAD

% Prob.ADObj = 1;
% madinitglobals;
% Result3 = tomRun('path',Prob,1);

% Approximate times:
% 1. Only f given:        0.25 s
% 2. Both f anf J given:  0.03 s
% 3. Only f with MAD:     0.11 s

% Adding a linear constraint.

% A = [0 0 0 0 0 0 0 0 0 1 1];
% b_U = 2;
% 
% Prob3 = mcpAssign('mcptest1_f', 'mcptest1_J', [], [], [], 'MCP 3', ones(11,1), A, [], b_U);
% Prob3.PriLevOpt = 1;
% 
% Result3 = tomRun('path',Prob3,1);