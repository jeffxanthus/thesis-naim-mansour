% Run some tests

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.0.0$
% Written Oct 12, 2000.    Last modified Mar 30, 2005.

%Result = tomRun('ucSolve','uc_prob',1,3);

%Result = tomRun('ucSolve',1,1,3);

%Result = tomRun('ucSolve',1,2,3);

%Result = tomRun('ucSolve',3,2,3);
%Result = tomRun('ucSolve',7,2,3);
%Result = tomRun('ucSolve',8,2,3);

%Result = tomRun('lpSimplex',1,2,3);

%Result = tomRun('lpSimplex',3,4,3);

%Result = tomRun('ucSolve');

%Result = tomRun('ucSolve',[]);
%Result = tomRun('ucSolve',[],[],[]);
%Result = tomRun('','uc_prob',1,[],1,3);
Result = tomRun('','con_prob',4,[],1,3);