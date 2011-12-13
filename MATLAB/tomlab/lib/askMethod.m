% askMethod calls SolverMethod to obtain a header and a menu,
% then asks the user to choose a method item from this menu.
%
% function method = askMethod(Solver,SolverAlg);
%
% Normally this method is the choice of how to solve some subproblem in
% the Solver, like the solution of an equation system, the solution of
% an overdetermined system or how to solve a quadratic programming subproblem
 
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2000-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.2.0$
% Written Nov 29, 1998.     Last modified Jun 6, 2008.

function method = askMethod(Solver,SolverAlg)

[mHeader,mMenu] = SolverMethod(Solver,SolverAlg);

if isempty(mHeader)
   method=0;
   fprintf('\n');
   fprintf('------------------------------------------------------------\n');
   fprintf('The solver %s  with SolverAlg = %d needs no method menu\n',...
            Solver, SolverAlg);
   fprintf('------------------------------------------------------------\n');
   fprintf('\n');
else
   method=-1+strmenu([Solver ': ' mHeader],mMenu);
end

% MODIFICATION LOG
%
% 981129  hkh  Written
% 080606  med  Cleaned xGUI