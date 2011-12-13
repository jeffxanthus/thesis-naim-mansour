% psfDemo.m
%
% Example showing the use of partially separable functions (psf)
%
% psf is only used by solver sTrustr
%
% The test example is #14 in con_prob, DAS2 from Hoch-Schittkowski

probFile = 'con_prob';              % Name of Init File
P        = 14;                      % Problem number in con_prob
Prob     = probInit(probFile, P);   % Define a problem structure

Prob.PriLevOpt = 3;

Result   = sTrustr(Prob);           % Call sTrustr
pause

% Notice that the number of function evaluations get #psf times larger
PrintResult(Result,2);

% Call conSolve as a comparison. This example is very easily solved
Result   = conSolve(Prob);
PrintResult(Result,2);
pause


% The only thing that triggers the use of the partial separability is
% the definition of the variable \VAR{Prob.PartSep.pSepFunc}.
% To solve the same problem, and avoiding the use of {\it psf},
% the following statements could be used:

probFile = 'con_prob';              % Name of Init File
P        = 14;                      % Problem number in con_prob
Prob     = probInit(probFile, P);   % Define a problem structure

Prob.PartSep.pSepFunc = 1;          % Redefining number of separable functions

% An alternative is to set the field as empty
%    Prob.PartSep= [];        

Prob.PriLevOpt = 3;

Result   = sTrustr(Prob);

PrintResult(Result,2);

% For this simple example there is no advantage in solving the problem
% using the psf facility. The additional chain of function calls just
% takes a lot of extra time.