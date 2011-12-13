% TOMLAB /CPLEX - LP/QP/MIP/MIQP/MIQQ Solver
% Version 7.7 (R7.7.0) 9-May-2011 
%
% Contents.m     This file
%
% startupcplex.m Using Matlab, go to TOMLAB /CPLEX directory, run
%                startupcplex.
%                Example:  cd \tomlab\cplex; startupcplex
%                The larger general TOMLAB toolbox can automatically detect
%                TOMLAB /CPLEX and initialize it.
%
%                Example:
%                cd d:\tomlab; startup
%                cd d:\tomlab\cplex; startupcplex
%
% install.m      Help in installing TOMLAB /CPLEX on PC systems
%
% install.unx    Help in installing TOMLAB /CPLEX on UNIX/Linux systems
%
% changes.m      Log of changes in TOMLAB /CPLEX
%
% ---------------------------
% Files in the main directory
% ---------------------------
%
% cplex.m       Main Matlab routine that calls CPLEX (non-Tomlab format)
% cplexTL.m     Converts the Tomlab Prob structure format to
%               a standard call to the cplex.m main CPLEX routine.
%               Called from mexMIP.m. Do NOT call mexMIP.m directly.
%
% cplexmex.dll  TOMLAB /CPLEX DLL calling the CPLEX solver. Extension may vary
%               depending on platform.   
%
% When calling CPLEX from Tomlab, use the driver routines:
%      tomRun('cplex', ...)
%
% -------------------------------------
% CPLEX Standard Distribution files
% -------------------------------------
% cplex*.dll   In tomlab\shared. Extension and exact name may vary depending on
%              platform.
%
% ------------------
% Callback routines:
% ------------------
% cpxcb_BARRIER.m    Barrier callback
% cpxcb_DISJCUT.m    Disjunctive cut callback
% cpxcb_DUAL.m       Dual simplex callback
% cpxcb_DUALCROSS.m  Dual crossover callback
% cpxcb_FLOWMIR.m    Mixed integer rounding cut callback
% cpxcb_FRACCUT.m    Gomory fractional cut callback
% cpxcb_INCUMBENT.m  MIP Incumbent callback
% cpxcb_MIP.m        MIP callback
% cpxcb_MIPPROBE.m   MIP probe or clique merging callback
% cpxcb_PRESOLVE.m   Presolve callback
% cpxcb_PRIM.m       Primal simplex callback
% cpxcb_PRIMCROSS.m  Primal crossover callback
% cpxcb_QPBARRIER.m  QP Barrier callback
% cpxcb_QPSIMPLEX.m  QP Simplex callback
% cpxcb_USERCUT.m    MIP User Cut callback
%
% The above files are copied and edited by the user for advanced callback
% handling.
%
% ---------
% Utilities
% ---------
% cpxPrint.m    Print the result of a CPLEX run
%
% cpx2cbinfo.m  Create global vector cpxCBInfo. Available at each callback
%               and after the run.
%               The description for each element is described in the Cplex
%               Users Guide and when doing help cpx2cbinfo.
% cpx2retvec.m  Create global vector cpxRetVec. Available after the run.
%               The description for each element is described in the Cplex
%               Users Guide and when doing help cpx2cbinfo.
%
% cpx2mat.m     Import an MPS file (and others) to MATLAB.
%
% cpx2matmex.dll MEX used by cpx2mat. Extension may vary depending on platform.
%
% ================
%
% ..............
%
% examples   A set of example files:
% ========
%
% cpxaircrew Solve an air-crew scheduling problem. Subfunctions used:
%            generateToDs:  Generates feasible ToDs (Tour of Duty)
%            sectordata:    Generates the data
%
% cpxbiptest Test of three large binary integer programs for different
%            values of the control parameters
% cpxiptest  Test of three large integer programs for different
%            values of the control parameters
% cpxtomtest1 Test solving LP, QP and MIP problems predefined in the
%            TOMLAB Init File (IF) format.
% cpxtomtest2 Solving simple LP, MIP problem using TOMLAB TQ format
%
% cpxKnaps   Test of three knapsack problems. Use different cut strategies.
%            Calls cplex.m directly.
% cpxKnapsTL The same routine as cpxKnaps, but using the TOMLAB (TQ)
%            input format and calling the Tomlab driver:
%                  Result = tomRun('CPLEX',Prob);
% cpxTest1   Generalized assignment problem from Wolsey 1998, 9.8.16, pp165.
%            Define the linear sos1 constraints explicitly.
%            No presolve, aggressive cut strategy.
% cpxTest2   Generalized assignment problem from Wolsey 1998, 9.8.16, pp165.
%            Define sos1 constraints, see cplex.m
%            No presolve, aggressive cut strategy.
% cpxTest3   Generalized assignment problem from Wolsey 1998, 9.8.16, pp165.
%            No presolve, no aggressive cut strategy.
% cpxTestConflict Demonstration of the TOMLAB /CPLEX Conflict Refinement
%                 feature
% cpxTestQP1 Simple quadratic programming (QP) problem, also solved as MIQP
% cpxTestQP2 Simple mixed-integer quadratic programming (MIQP) problem
