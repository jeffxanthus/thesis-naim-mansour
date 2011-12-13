echo('on');

% ----------------------------------------------------------------------------
% Run a constrained problem, without defining or using an Init File and a GUI.
% ----------------------------------------------------------------------------
%
% Then the TOMLAB Quick format should be used, to setup a Prob problem structure
%
% There are a set of different ...Assign routines for different problem types
%
% Best way is to use the routine conAssign to define the Prob structure for
% constrained nonlinear problems.
% -------------------------------------------------------------------------

Prob   = conAssign('con_f','con_g','con_H',[], [-100 -100],[],'con test',...
            [.1 2],[],[],[],[],[],'con_c','con_dc',[],[], zeros(2,1));

Result = tomRun('conSolve',Prob,2);  % Solve the problem, and print result

% -------------------------------------------------------------------------
% The con_prob Init File have many nonlinear problems defined. 
% If to run one of these problems from the command line, 
% the following is the way, if to be able to change some input parameters
% -------------------------------------------------------------------------

Prob = probInit('con_prob',1);          % Select problem 1
Prob.optParam.IterPrint = 1;            % Output every iteration
Result = tomRun('conSolve',Prob,2);  % Solve the problem, and print result

% ----------------------------------------------
% Otherwise one can just run a problem directly:
% ----------------------------------------------

Result = tomRun('nlpSolve','con_prob',2);
echo('off');
