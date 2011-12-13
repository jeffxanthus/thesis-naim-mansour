% askAlgorithm calls SolverAlgorithm to obtain a header and a menu,
% then asks the user to choose an algorithm item from this menu.
%
% function SolverAlg = askAlgorithm(Solver);

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2000-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.2.0$
% Written Sept 25, 2000. Last modified Jun 6, 2008.

function SolverAlg = askAlgorithm(Solver)

[mHeader,mMenu] = SolverAlgorithm(Solver);

if isempty(mHeader)
   SolverAlg=0;
   fprintf('\n\n');
   fprintf('------------------------------------------------------------\n');
   fprintf('The solver %s needs no choice of algorithm\n',...
            Solver);
   fprintf('------------------------------------------------------------\n');
   fprintf('\n\n');
else
   SolverAlg=-1+strmenu([Solver ': ' mHeader],mMenu);
end

% MODIFICATION LOG
%
% 000925  hkh  Written
% 080606  med  Cleaned xGUI