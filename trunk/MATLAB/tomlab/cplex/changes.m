% LOG OF CHANGES IN TOMLAB /CPLEX
%
% December 2007:
%
% TOMLAB /CPLEX 11.0.0
% 
% TOMLAB /CPLEX - CPLEX 11 now embedded in the official release.
%               - New tuning tool that helps the user minimize the execution time in production.
%               - Parallel mode possible for faster solution of mixed-integer models.
%               - Possible to collect many solutions in a pool.
%               - A wider variety of nonconvex quadratic constraints are now automatically handled.
%               - No longer printing error message when empty linear constraint matrix.
%               - CPLEX no longer crashing when providing xIP for MILP/MIQP problems.
%               - Using embedded license by default.
%               - Branching priority added to input variables.
%               - Possible to give branching directions for individual variables.
%               - MIP Incumbent and User Cut callbacks
%
% February 2007:
%
%   Fixed problem with QP/MIQP/MIQQ with all zeros quadratic matrix
%
% December 2006:
%
%   String buffer bug fixed.
%
% September 2006:
%
%   Screen printing improved to avoid delayed output 
%
% May 2006:
%
%   Unsymmetric QP handled more smoothly. 
%
% February 2006:
%
% TOMLAB /CPLEX 10.0.0
%
% Based on ILOG CPLEX 10.0.0
%
% New features:  Conflict Refinement (replaces IIS)
%                Indicator Constraints
%                More useful callbacks (see cpxcb_*.m)
%
% -----------------------------------------------------------------------
%
% TOMLAB /CPLEX 9.0.2
%
% Based on ILOG CPLEX 9.0.2
%
% New features: IIS (Irreducible Infeasible Sets)
%               SA  (Sensitivity Analysis)
%
% -----------------------------------------------------------------------
%
% TOMLAB /CPLEX 9.0
%
% 040204 TOMLAB /CPLEX 9.0 released, based on ILOG CPLEX 9.0.
%
%        Support for Quadratic Constraints.
%
% 040416 Network solver added.
%
% -----------------------------------------------------------------------
%
% TOMLAB /CPLEX 8.1.0
%
% 031024 TOMLAB /CPLEX released, now based on ILOG CPLEX 8.1.
%
%        Several new features added:
%
%          Possible to save the CPLEX problem in 6 different formats
%          Improved printing to the Matlab command window
%          LogFile parameter added
%
% -----------------------------------------------------------------------
%
% TOMLAB /CPLEX 8.0.1
%
% 021030 Added SAV write feature if PriLev >= 1000
%
% 020930 RELEASE OF TOMLAB /CPLEX v1.0, based on Release 8.0
%
% 020921 Major revision to make CPLEX examples similar to Xpress examples
%
% 030801 Fixed ExitFlag handling in cplexTL to conform with Tomlab standard
%
help changes