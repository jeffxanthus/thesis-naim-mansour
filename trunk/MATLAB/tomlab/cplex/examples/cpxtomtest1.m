% Test of TOMLAB /CPLEX call using TOMLAB input format
%
% Test of problems predefined in TOMLAB.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written July 28, 2001.     Last modified Aug 6, 2009.

% First solve a simple LP problem

Prob = probInit('lp_prob', 1);
Result = tomRun('cplex',Prob,1);

disp('Press ENTER to continue')
pause

Prob = probInit('qp_prob', 1);
Result = tomRun('cplex',Prob,1);

disp('Press ENTER to continue')
pause

Prob = probInit('miqp_prob', 1);
Result = tomRun('cplex',Prob,1);

disp('Press ENTER to continue')
pause

Prob = probInit('mip_prob', 14);
Prob.MIP.cpxControl.PREIND = 1;
Result = tomRun('cplex',Prob,14);

disp('Press ENTER to continue')
pause

Prob = probInit('mip_prob', 15);
Result = tomRun('cplex',Prob,15);

% MODIFICATION LOG:
%
% 020922 hkh Revised for CPLEX
% 050113 med Removed trial print-outs
% 050209 med Changed name from tomtest1 to cpxtomtest1
% 050829 med Changed to miqp prob call
% 090805 med Updated tomRun calls