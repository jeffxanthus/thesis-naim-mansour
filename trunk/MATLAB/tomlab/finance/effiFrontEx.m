% effiFrontEx.m:
%
% Example for Mean Variance Efficient Frontier using TOMLAB.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2005 by Tomlab Optimization Inc., $Release: 4.6.0$
% Written Dec 18, 2004.   Last modified Dec 19, 2004.

clear all

load effiFrontExMat
% Set UseFrontcon to 0, 1 or 2
% 0 = Only run effiFront, 1 = run effiFront + opt tbx frontcon, compare results
% 2 = Run effiFront + TOMLAB frontcon, compare results
UseFrontcon = 0;
if UseFrontcon == 1
    rmpath('c:\tomlab\optim');
    % C = 0.5*(C + C');  % Make C symmetric
    tic
    [PortRisk, PortReturn, PortWts] = frontcon(Er, C, 24, [], AssetBounds, groups_f, groups_f_bounds);
    toc
    % If to plot, uncomment this
    % frontcon(Er, C, 24, [], AssetBounds, groups_f, groups_f_bounds)
end

if UseFrontcon == 2
    addpath('c:\tomlab\optim');
    %C = 0.5*(C + C');  % Make C symmetric
    tic
    [PortRisk, PortReturn, PortWts] = frontcon(Er, C, 24, [], AssetBounds, groups_f, groups_f_bounds);
    toc
    % If to plot, uncomment this
    % frontcon(Er, C, 24, [], AssetBounds, groups_f, groups_f_bounds)
end

if UseFrontcon == 0
    PLOT = 0;             % 0 = no plot, 1 = plot efficient frontier
    Prob.PriLevOpt = 0;   % Print level in tomRun, 0 = silent, 1 = some printing

    % Could set alternative solvers in effiFront
    if 0
        Prob.SolverLP = 'qpopt';   % Set explicit solver for LPs
        Prob.SolverQP = 'qpopt';   % Set explicit solver for QPs
    elseif 0
        Prob.SolverLP = 'xpress-mp';
        Prob.SolverQP = 'xpress-mp';
    elseif 0
        Prob.SolverLP = 'cplex';
        Prob.SolverQP = 'cplex';
    elseif 0
        Prob.SolverLP = 'sqopt';
        Prob.SolverQP = 'sqopt';
    end

    [pRisk, pRet, pWts] = effiFront(Er, C, 24, [], PLOT, ...
        AssetBounds(1,:), AssetBounds(2,:), ...
        groups_f, groups_f_bounds(:,1), groups_f_bounds(:,2), Prob);

    if UseFrontcon > 0
        fprintf('\n');
        fprintf('Sum of difference in Risk values %f\n', sum(pRisk - PortRisk));
        fprintf('Sum of difference in Return values %f\n', sum(pRet - PortReturn));
        fprintf('\nNegative values above means Tomlab found better optima\n\n');
        fprintf('Sum of difference in each set of portfolio weights\n');
        fprintf('\n');
        sum(pWts - PortWts')
    end
end

% MODIFICATION LOG:
%
% 041221  hkh  Written