% nllsQG is a small example problem for defining and solving
% nonlinear least squares using the TOMLAB format.
Name='Gisela';

t = [0.25; 0.5; 0.75; 1; 1.5; 2; 3; 4; 6; 8; 12; 24; 32; 48; 54; 72; 80;...
     96; 121; 144; 168; 192; 216; 246; 276; 324; 348; 386];
y = [30.5; 44; 43; 41.5; 38.6; 38.6; 39; 41; 37; 37; 24; 32; 29; 23; 21;...
     19; 17; 14; 9.5; 8.5; 7; 6; 6; 4.5; 3.6; 3; 2.2; 1.6];
   
x_0 = [6.8729, 0.0108, 0.1248]';

% See help clsAssign for more information.
Prob = clsAssign('nllsQG_r', 'nllsQG_J', [], [], [], Name, x_0, ...
                 y, t);

% Parameter which is passed to r and J routines.             
Prob.uP = 5;

Result = tomRun('clsSolve', Prob, 1);
%Result = tomRun('nlssol', Prob, 1);

% Any nonlinear solver can be used. TOMLAB automatically
% uses gateway routines for problem mapping.

%Result = tomRun('filterSQP', Prob, 1);
%Result = tomRun('knitro', Prob, 1);
%Result = tomRun('conopt', Prob, 1);
%Result = tomRun('snopt', Prob, 1);
%Result = tomRun('npsol', Prob, 1);
%Result = tomRun('minos', Prob, 1);
%Result = tomRun('oqnlp', Prob, 1);