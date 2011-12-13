% Test of TOMLAB /CPLEX solver tuning capabilities.
%
% function [cpxControl1, cpxControl2] = cpxSolverTuning
%
% Exemplifies the use of the solver tuning features in TOMLAB /CPLEX and
% the parameters associated with this, TUNINGDISPLAY, TUNINGREPEAT,
% TUNINGTILIM and TUNINGMEASURE.

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc., $Release: 11.2.0$
% Written Nov 26, 2008.   Last modified Nov 26, 2008.

function [cpxControl1, cpxControl2] = cpxSolverTuning

% Create a test problem from the TOMLAB sets
Prob = probInit('mip_prob', 18);

% Choose to tune, but not to solve the problem
% A setting of 1 also solves the problem
Prob.CPLEX.Tune = 2;

% Display standard report plus parameter settings being tried
Prob.MIP.cpxControl.TUNINGDISPLAY = 2;

% Repeat the tuning 2 times for model
Prob.MIP.cpxControl.TUNINGREPEAT = 2;

% Tune time limit of 5 seconds (for each test)
% Use TILIM if limiting the entire run
Prob.MIP.cpxControl.TUNINGTILIM = 5;

% Set mean average of time to compare different parameter sets
Prob.MIP.cpxControl.TUNINGMEASURE = 1; % Default
Prob.PriLevOpt = 1;
Result = tomRun('cplex', Prob, 0);
cpxControl1 = Result.MIP.cpxControl;

% Minmax approach to compare the time of different parameter sets
Prob.MIP.cpxControl.TUNINGMEASURE = 2;
Prob.PriLevOpt = 1;
Result = tomRun('cplex', Prob, 0);
cpxControl2 = Result.MIP.cpxControl;

fields = {'TUNINGDISPLAY', 'TUNINGREPEAT', ...
    'TUNINGTILIM'};
cpxControl1 = rmfield(cpxControl1,fields);
fields = {fields{:},'TUNINGMEASURE'};
cpxControl2 = rmfield(cpxControl2,fields);

disp(' ');
disp(sprintf('Option 1: A total of %d parameters were modified.', ...
    length(fieldnames(cpxControl1))));
disp(' ');

disp(sprintf('Option 2: A total of %d parameters were modified.', ...
    length(fieldnames(cpxControl2))));

% MODIFICATION LOG:
%
% 081126 med Written