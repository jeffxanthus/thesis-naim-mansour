% SolverMethod returns a header and a menu for a particular solver 
%
% The different items in the menu is the different choices of method
% for some subtask in the given solver, most often how to solve
% a linear equation system, an overdetermined linear system, or
% a quadratic programming subproblem
%
% function [Header, Menu] = SolverMethod(Solver,SolverAlg)
% 
% INPUT: 
% Solver      Name of solver
% SolverAlg         A number selecting one of the algorithms for the Solver.
%
% OUTPUT: 
% Header      The header
% Menu        A string matrix with the different menu options

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., $Release: 5.2.0$
% Written Oct 28, 1998. Last modified Feb 12, 2006.

function [Header, Menu] = SolverMethod(Solver,SolverAlg)

if nargin < 2
   SolverAlg=[];
end

if isempty(SolverAlg), SolverAlg=0; end

Header=[]; Menu=[];

switch lower(deblank(Solver))
   case 'ucsolve'
      if any(SolverAlg==[0 1 2 4])
         Menu = str2mat( ...
            'Singular value decomposition' ...
           ,'LU-decomposition' ...
           ,'LU with pivoting' ...
           ,'Left matrix divide (QR)' ...
           ,'Left matrix divide' ...
           ,'Explicit inverse (inv)' ...
         );
         Header='Linear Equation Method';
      elseif any(SolverAlg==[6 7 8])

         % If c-g method
         Menu = str2mat( ...
            'No reset' ...
           ,'Reset after n iterations' ...
         );
         Header = 'Reset technique';
      end
   %case {'qpsolve'}
   %case {'consolve'}
   %case {'nlpsolve'}
   %case {'strustr'}
   case {'clssolve'}
      Header='Linear Equation Method';
      %Header= 'GN Step Solution Method';
      Menu = str2mat( ...
         'QR-decomposition' ...
        ,'Singular value decomposition' ...
        ,'Matlab built-in inverse' ...
        ,'Matlab pseudoinverse' ...
      );

   case {'mipsolve'}
      Header = 'Variable selection';
      Menu = str2mat(...
         'Simplex, Minimum cost rule' ...
        ,'Simplex, Blands rule' ...
        ,'Simplex, MC rule (Dantzig)' ...
      );
   case {'cutplane'}
      Header = 'Variable selection';
      Menu = str2mat(...
         'Simplex, Minimum cost rule' ...
        ,'Simplex, Blands rule' ...
        ,'Simplex, MC rule (Dantzig)' ...
      );
   %case {'lpsimplex'}
   %case {'glbsolve'}
   %case {'glbfast'}
   %case {'glcsolve'}
   %case {'glcfast'}
   %case {'glccluster'}
   %case {'rbfsolve'}
   %case {'arbfmip'}
   %case {'ego'}
   %case {'dualsolve'}
   %case {'infxolve'}
   %case {'lpopt','qpopt','minos''npsol','nlssol','snopt','sqopt'}
   %case {'fmincon'}
   %case {'fminsearch'}
   %case {'fminunc'}
   %case {'linprog'}
   %case {'lsqcurvefit'}
   %case {'lsqlin'}
   %case {'lsqnonlin'}
   %case {'lsqnonneg'}
   %case {'quadprog'}
   %case {'constr'}
   %case {'fmins'}
   %case {'fminu'}
   %case {'lp'}
   %case {'qp'}
   %case {'leastsq'}
   %case {'sqopt'}
   %case {'qld'}
   %case {'lssol'}

   %case {'xpress-mp'}
   %case {'bqpd'}
   %case {'minlpbb'}
   %case {'miqpbb'}
   %case {'filtersqp'}
   %case {'cplex'}
   %case {'pensdp'}
   %case {'penbmi'}
   %case {'pdco'}
   %case {'pdsco'}
end

% Use the following, or???
%,'Affine Scaling Karmarkar' ...
%,'lpsimp2, Minimum cost rule' ...
%,'lpsimp2, Blands rule' ...
%,'lpsimp2, MC rule (Dantzig)' ...


% MODIFICATION LOG
%
% 981129   hkh  Routine written
