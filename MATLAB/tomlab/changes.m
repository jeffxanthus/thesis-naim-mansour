% LOG OF MAJOR CHANGES IN TOMLAB
%
% Please refer to: http://tomopt.com/tomlab/company/news.php
% for information about the latest release.
%
% 091030 TOMLAB 7.3 released
% 090818 TOMLAB 7.2 released
% 090325 TOMLAB 7.1 released
% 081117 TOMLAB 7.0 released
%
% TOMLAB now features a complete modeling platform with source
% transformation: TomSym
%
% 080606 TOMLAB 6.1 released
%
% TOMLAB /BASE - Version 6.1 released.
%        - glcDirect printing updated.
% 		 - glcCluster using less memory for problems with integer variables.
% 		 - mipSolve using subsolver in the following order (determined by license): MINOS, BQPD, XA, MILPSOLVE.
% 		 - Substantial updated to Linux distributions (library dependencies and compilers).
% 		 - glcCluster updated to use maxLocalTry.
% 		 - Problems when running quadratically constrained problems of varying sizes fixed.
% 		 - nlpSolve updated to handle lp subproblems with an objective.
% 		 - NumDiff and ConsDiff automtically set to correct values when empty.
% 		 - Second order information not computed for problems without analytical derivatives.
% 		 - nlpSolve now has 500 iterations as default.
% 		 - Optimization toolbox interface updated.
% 		 - Now poossible to use MAD with multiMin.
% 		 - linRatSolve automatically purging binary and integer variables.
% 		 - Optimization toolbox interface updated.
% 		 - miqqAssign updated to support problems without quadratic constraints.
% 		 - bincont2lin now taking lower bounds into account.
%
% TOMLAB /CPLEX - MAC OS X (Intel) now supported.
% 		  - Advanced presolve functionality enabled.
% 		  - Help updated.
% 		  - Problem with incorrect objective reported fixed.
% 		  - b_L and b_U verified to be dense.
% 		  - Change obsolete parameter MIPSTART into ADVIND.
% 		  - Screen printing disabled when using more than 1 thread.
%
% TOMLAB /MAD - The following functions added:
% 		     FFT, FFTN, IFFT, FFTSHIFT, IFFTSHIFT.
% 		- New website for MAD: http://matlabad.com/
% 		- TRANSPOSE updated for improved speed.
% 		- IMAG and REAL added.
%             - General updates to derivvec class.
% 		- MIN updated.
%
% TOMLAB /KNITRO - Version 5.2 embedded.
% 	         - New options BAR_PENALTYCONS and BAR_PENALTYRULE
% 	         - New options for multi-start feature: MS_NUMTOSAVE, MS_SAVETOL and MS_TERMINATE
% 	         - New options for file printing: OUTDIR and OUTAPPEND
% 		   - New change for FEASIBLE and FEASMODETOL to BAR_FEASIBLE  and BAR_FEASMODETOL.
%
% TOMLAB /MINLP
%       - Improved robustness when error conditions occur in user code
%
% TOMLAB /PENOPT - PENSDP and PENBMI added for Linux 64-bit.
%
% TOMLAB /SOL - Defaults updated for LPOPT, QPOPT, SNOPT7, SNOPT, SQOPT, SQOPT7.
%
% TOMLAB /GENO - Fixed bug when x_0 present.
%
% TOMLAB /AMPL - Writing solution file regardless of ExitFlag.
%
% TOMLAB /PROPT - Released (Lead Developer is Per Rutquist - Matlab Programming Contest champion).
% 		  - A complete framework for optimal control (ODE and DAE) and parameter estimation.
% 	      - Over 100 fully functional test cases.
% 		  - See the dedicated webpage here: http://tomdyn.com/
%         - 6 years of optimal control experience has proven the software unbeatable in comparison to competing modules.
%
% 071228 TOMLAB 6.0 released
%
% TOMLAB /BASE - Version 6.0 released.
%              - All short-circuit and/or statements reversed to standard | and &.
%              - It is now possible to define the problem with the following functions:
%                  Nested function
%                  Anonymous function
%                  Subfunctions defined within main file
%                  Regular m-files referenced by function handle or string
%              - glcDirect error reporting updated when using function handles.
%              - Tfmin updated and added to manual.
%              - Several Base Module solvers updated to handle non-symmetric matrices.
%              - glcCluster warm start updated.
%              - Tfzero help updated (not possible to use gradient) and added to manual.
% 	            - glcCluster now using Prob.x_0 and Prob.X0 as initial starting points for local searches.
% 	            - multiMin - new default. Now using min(3000,30*n) local searches.
% 	            - ucDemo updated to TOMLAB format.
% 	            - tomGui and other GUIs removed from distribution.
% 	            - Major updates to the general TOMLAB manual.
% 	            - Manual updates for MAD, CPLEX, KNITRO, MINLP, SOL.
%
% TOMLAB /SOL - Updated SOL manuals with regards to basis files.
%
% TOMLAB /KNITRO - Problem removed when no nonlinear constraints are present.
%                - Function evaluation report corrected.
%
% TOMLAB /AMPL - Updated to handle LCP, QCP and MCP problems.
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
% TOMLAB /MAD - Updated interp routines to work for most future Matlab releases.
%               AD for second order constraints working.
%
% TOMLAB /PENOPT - Now works with Matlab 7.5.
%
% TOMLAB /CONOPT - Available for Windows 64-bit.
%                - Printing removed when empty linear constraints are included.
%
% TOMLAB /OQNLP - LSGRG2 printing updated.
%
% TOMLAB /CGO - Bug when using glcCluster as subsolver removed.
% 	    - rbfSolve and arbfmip (not released yet) updated for problems with stationary sample points.
% 	    - Three new experimental designs (-997, -998, -999 (revised)). Strategy -999 is no longer called the Gutmann strategy.
% 	    - Updated when using multiMin as global solver.
% 	    - New REPLACE option for rbfSolve and arbfmip.
%      - Max Lipschitz constant LipMx estimated using the initial set of points X and updated in every iteration.
%      - New defaults for box-constrained problems (REPLACE = 5, SCALE = 0).
% 	    - AddSurfMin new parameter for checking additional interior local minima.
%      - New parameter when using multiMin.
% 	    - RandState modified for EGO when no warm start.
%
% 070907 TOMLAB 5.9 released
%
% TOMLAB /BASE - Version 5.9 released.
%                Bug in PrintAssign fixed.
%                Several updates to lsqcurvefit.
%                conSolve updated gradient handling.
%                Compatible with Matlab 7.5 (R2007b).
%                clsSolve updated for Broydens method - Jacobian set correctly initially.
%                Bug in qpSolve removed.
%                ucSolve printing updated.
%                glcDirect updated to take 1 or 2 inputs.
%                multiMin updated to enable warm start of solution process.
%
% TOMLAB /SNOPT - SNOPT 7.2-5 (May 2007)
%
% TOMLAB /KNITRO - Memory bug fixed.
%                - Improved Hessian handling for large problems.
%
% TOMLAB /CPLEX - New tool for adding piece-wise linear functions and constraints to an existing model.
%               - User cuts can now be used in non-parallel mode.
%               - Incumbent callback available.
%
% TOMLAB /XA - Help updated.
%
%
% 070625 TOMLAB 5.8 released
%
%
% TOMLAB /BASE  - Version 5.8 released.
%               - milpsolve MaxCPU parameter added.
%               - Updates implemented for MAC OS X (Intel).
%               - Warm start code updated for most solvers.
%               - License not found print-out from MAD removed.
%               - Possible to generate a print file from TLSQR.
%               - estBestHessian updated for use with MAD.
%               - License problems for Windows 2000 and Windows Vista corrected.
%
%
% TOMLAB /CGO   - Safer memory handling in tomsol CGO support library
%
%
% TOMLAB /CPLEX - Problem with objective for linear problems resolved.
%               - Conflict refiner can now generate a print file.
%               - Callback added: cpxcb_INCUMBENT called for all new integer solutions.
%               - Review of all return and error codes.
%               - Now based on CPLEX 10.2
%
%
% TOMLAB /XA    - Minor incompatibility problem fixed with other toolboxes.
%
%
% TOMLAB /SOL   - NPSOL return code and text modified.
%               - LPOPT and QPOPT updated for warm start.
%               - NPSOL updated for automatic differentiation of nonlinear constraints.
%
%
% TOMLAB /KNITRO - Updates to manual.
%
%
% TOMLAB /CONOPT - Notation updated to match standard.
%                - Updated to handle empty constraint Jacobian.
%
%
% TOMLAB /OQNLP - New version with many bug fixes.
%
%
% TOMLAB /GP    - Memory problem corrected.
%
%
% TOMLAB /MAD   - Bug in Matlab prevents certain object oriented code from working properly
%                 (avoid Matlab 7.4 with MAD)
%               - Added function: tanh, sinh, cosh, coth, atanh, acoth and asinh.
%               - userfunc.m added as template for user functions.
%               - interp1 and interp2 compatible with Matlab 7.4.
%
% 070216 TOMLAB 5.7 released
%
% TOMLAB /BASE   - Version 5.7 released.
%                - MATLAB R2007A compatible
%                - ADMAT removed from distribution
%
% TOMLAB /CPLEX  - Improved stability for multiple threads in conjunction with
%                  screen output.
%
%                - Correctly handles QP/MIQP/MIQQ problems with zero quadratic
%                  objectives
%
%
% 061205 TOMLAB v5.6 released
%
% http://tomlab.biz has now been changed to http://tomopt.com. Please update your bookmarks!
%
% TOMLAB /BASE   - Version 5.6 released.
%                - Separate 64-bit versions developed for Matlab 7.3 (R2006b) and later.
%                - Several updates the TOMLAB test suite.
%                - Writing to MAT-files from global solvers (glb/glcDirect, etc) safer, avoiding losing results after long runs.
%                - Minor help updates to m-files.
%                - qpSolve updated to work with nonlinear subsolvers.
%
% TOMLAB /MAD    - interp1, interp2 and max updated.
%                - not added to distribution.
%                - cat added - possible to do horizontal and vertical concatenation of MAD objects.
%                - general bug fixes and updates.
%
% TOMLAB /CPLEX  - Bug related to string buffer fixed.
%
% TOMLAB /CGO    - Writing to MAT-files safer, avoiding losing results after long runs.
%
% TOMLAB /KNITRO - MAXTIME parameter corrected.
%                - Updated to handle problems with one dense nonlinear constraint.
%
% TOMLAB /SOL    - SNOPT and SNOPT7 no longer making callbacks for problems with pure linear objectives.
%                - MINOS updated to handle nonlinear problems with linear objectives.
%                - Differences between NPSOL for 32-bit and 64-bit Windows resolved.
%
% TOMLAB /MSNLP  - MSNLP and LSGRG2 Now available for 64-bit Linux.
%                - Stability improved on all platforms.
%
%
% 060905 TOMLAB v5.5 released
%
% TOMLAB /BASE - Version 5.5 released.
%                Issues with number of inputs for glcDirect and glcCluster resolved.
%                chs_prob updated.
%                Updates to InstallShield installer for Windows (32/64 bit).
%                Several manual updates.
%                Functions for modifying static parts of existing problems added (modify_* and more).
%                  See section 4.5 in the main TOMLAB manual.
%                LDO and LDODEMO removed from distribution (available for download from the manuals page).
%                License problem with Matlab R12 resolved.
%                lpSimplex no longer printing by default.
%                TOMLAB version of LSQCURVEFIT bug fixed.
%                clsSolve - Tests added to avoid running problems with nonlinear constraints.
%                         - Revision of Broydens method (line search, defaults)
%                         - No longer crash when Jacobian is all zeros.
%                goalSolve - Updates for handling of nonlinear constraints.
%                          - Handling of Lagrange multipliers updated.
%                Improvements maybe to optimization toolbox plug-in FGOALATTAIN.
%                Tolerances updated for printing.
%                PDSCO now working with problems without linear constraints.
%                Output improved when problems arise.
%                System-wide updates for defaults in assign routines.
%                infSolve and infLinSolve updated for lower bound setting, now possible to get an objective lower than 0.
%                Updates to defaults for line searches used.
%                Consistent handling of convergence to user defined target value.
%                Unnecessary isnan checks removed from all solvers.
%                LINPROG plug-in print outs removed.
%                multiMin adjusted to also handle Matlab 6 and least squares problems.
%                Additional safety added to all assign routines.
%                nlresid and nlfunc no longer calling nlp_c for problems with no nonlinear constraints.
%
%   TOMLAB /GENO - A package for genetic programming release.
%
%   TOMLAB /CGO  - Pre-defined test problems without IntVars can now be run with rbfSolve.
%                - Tests if bounds are finite.
%
%   TOMLAB /SOL  - Manuals updated with return codes for SQOPT.
%                   LSSOL/LPSOL/QPOPT/SQOPT updated for problems with no linear constraints.
%                   Screen printing buffers corrected for SNOPT7 and SQOPT7.
%                   SQOPT now accepting problem with empty linear part.
%
%   TOMLAB /MINLP - f_Low now used for the solvers.
%                   filterSQP - Lagrange multipliers now calculated correctly.
%                   Help updated regarding dense and sparse solver versions.
%
%   TOMLAB /MAD   - Several general updates and improvements to the system.
%                   interp1 and interp2 added to functionality.
%
% 060814 TOMLAB v5.4 released
%
%   TOMLAB /BASE  - Help for assign routines updated.
%                 - ucSolve help updated.
%                 - New solver multiMin that obtains multiple local minima from either a given or random set of starting points. A local sub-solver is required.
%                 - New simpler licensing system - substantial speed-ups for demo users.
%                 - TSP (Travelling salesman problems) integrated in tomlab\testprob.
%                 - sqr2 removed from distribution.
%                 - Updated optimization toolbox interfaces.
%                 - 21 Test problems added from GENO manual.
%                 - clsSolve now using TLSQR as subsolver for large-scale NLLS problems.
%                 - PDCO has been improved.
%                 - mipSolve convergence improvements.
%
%   TOMLAB /CPLEX - Screen printing fixed for CPU intensive problems (bug in Matlab prevented prints from being flushed).
%
%   TOMLAB /KNITRO - Screen printing fixed for CPU intensive problems.
%
%   TOMLAB /CONOPT - Screen printing fixed for CPU intensive problems.
%
%   TOMLAB /XPRESS - Screen printing fixed for CPU intensive problems.
%
%   TOMLAB /MINLP - Updated help.
%
% 060701 TOMLAB v5.3 released
%
%   TOMLAB - Version 5.3 released.
%        - No lower limit on integer solution set for glcDirect anymore.
%        - BuildMPEC added to the distribution. Support with the KNITRO solver now possible.
%        - glcSolve output improved.
%        - All pattern estimations updated (estConsPattern, HessPattern, d2cPattern, JacPattern).
%        - New quick guide with examples for linear, quadratic and nonlinear complementarity problems.
%        - Three new assign routines for complementarity problems.
%
%   TOMLAB /KNITRO  - Version 5.0 released
%                 - New multistart features.
%                 - Possibility to solve MPEC problems.
%                 - See version 5.2.1 for more information.
%
%   TOMLAB /CGO - Algorithmic convergence greatly improved.
%               Printing improved during optimization process.
%               Rescue proceduce implemented.
%               Bug fixed that could cause a crash.
%
%   TOMLAB /MAD - Complete system review.
%               Contents.m files added to all folders.
%               Facilities for differentiating through black box code.
%               High-level interfaces allowing easy use of MAD with Matlab's ODE solvers and many Optimization Toolbox solvers.
%               Users guide updated with more examples.
%
% 060609 TOMLAB v5.2.1 released
%
%     * TOMLAB /KNITRO 5.0 released.
%
%       Now based on KNITRO 5.0.2. New features/changes:
%
%          Improvements in efficiency and robustness across all 3 algorithms.
%          The most significant improvements were made in the "Active" algorithm.
%
%          A crossover procedure for "cleaning up" the interior-point solution
%          by switching to the active-set at the end has been added.
%
%          Mathematical Programs with Equilibrium Constraints (MPECs) can now be solved.
%
%          Basic heuristics for multi-start global optimization added.
%
%          New BARRULE setting (6/QUALITY) added.
%
% 060504 TOMLAB v5.2 released.
%
%        - Several manuals are now available in HTML as well as PDF.
%        - A new manual covering all the test problems included with TOMLAB released.
%        - All assign routines now checking for crossover bounds.
%        - Updated error checking for clsSolve.
%        - FMINUNC interface updated to handle problems with only an objective function.
%        - miqqAssign updated to set c_L and c_U correctly.
%        - ls_H updated to handle nonlinear least squares problems more efficiently.
%        - glcDirect/glcFast - input checking now fixed.
%        - glcCluster - major updates to algorithm.
%        - infSolve - now possible to include integer variables.
%
%   TOMLAB /NPSOL - No longer estimating ConsPattern for large scale problems.
%
%   TOMLAB /CGO - Changed glcFast to glcDirect as default, added glcDirect setup.
%               - New inputs nSample and eps_sn.
%               - Use new function expDesign for initial experimental design.
%
%   TOMLAB /SNOPT - SNOPT 7 added to warm start functionality.
%
%   TOMLAB /CPLEX - Unsymmetric QP problem handled more smoothly.
%
% 060206 TOMLAB v5.1 released.
%
%     * System wide manual updates.
%     * Routine lls2qp added to distribution. Now possible to convert a problem with an LLS objective to a problem with a QP objective.
%     * Now possible to solve mixed-integer linear least squares problem (test problem added to the Quickguide).
%     * expSolve (solver for fitting of sums of positive exponentials) separated into expAssign and expSolve.
%     * Possible to call slsSolve, L1Solve, infSolve with tomRun and get normal printing.
%
%
%   TOMLAB /CGO:
%
%     * A new version of EGO released - major updates to the entire algorithm stucture - the solver has proven more stable with the new c-libraries.
%     * DACEProb is no longer a copy of main Prob.
%     * Transformation for EGO callback function ego_f added. Can be controlled with EITRANSFORM.
%
%
%   TOMLAB /NPSOL:
%
%     * NLSSOL no longer doing isfield check in callback.
%
%
%   TOMLAB /SOL:
%
%     * Verify level default set to -1 for all solver. No derivative checks done.
%     * Informs values revised for all solvers. New manual and m-files available.
%
%
%   TOMLAB /MINLP:
%
%     * Branch strategies for miqpBB updated. Negative variables were causing problems before.
%
%
%   TOMLAB /CPLEX - Version 10 released.
%
%     * A new more advanced conflict refiner has replaced the IIS features.
%     * Possible to give logical (indicator) constraints (i.e. a binary variable controls the constraint).
%     * The following parameters have been added/removed/modified:
%
%           COVERS - Now possible to do very aggressive cover cuts.
%           PRESLVND - One more option added for node presolve.
%           SYMMERTY - Several new options to generate these cuts.
%           FEASOPTMODE - Added with many 5 possible settings.
%           REPEATPRESOLVE - Controls how to re-apply presolve for MIP models.
%           MEMORYEMPHASIS - Memory may be conserved using this parameter.
%           NUMERICALEMPHASIS - Caution parameter.
%           POLISHTIME - Time spent polishing MIP solution.
%           EPRELAX - Control for FEASOPTMODE.
%
%           BAROOC, FINALFACTOR, PRECOMPRESS, COLGROWTH, NZGROWTH, QPNZGROWTH, ROWGROWTH,
%           REVERSEIND, XXXIND and BASINTERVAL removed.
%
%
% 051215 TOMLAB v5.0 released:
%
%         - estd2cPattern updated (Prob.N was used instead of Prob.mNonLin)
%         - portfolioAbsCon updated.
%         - Updated error messages for installation problem.
%         - binbin2lin added to the distribution - automatically converts
%           problems with binary products. One can then model with an extra
%           variable right away without re-writing the problem.
%         - bincont2lin added to the distribution - converts problems with
%           binary integer/continuous products.
%         - New solvers added:
%
%           INFLINSOLVE - solves linear and mixed integer linear minimax problems.
%           LINRATSOLVE - solver linear programming problem with a ratio for the objective.
%
%           See the TOMLAB manual for more details.
%
%         - Two new problem types added to the TOMLAB Quickguide.
%         - New MaxIter, MaxFunc defaults for glbSolve, glbFast, glcSolve
%           and glcFast (MaxIter = max(5000,n*1000), MaxFunc = max(10000,n*2000)).
%         - Assign routines updated with more checks on inputs.
%         - mipSolve, cutplane and minlpSolve updated.
%
%         - Model library for linear and mixed-integer programming released.
%           Detailed test cases for the following areas:
%
%             - air transportation
%             - finance and economics
%             - ground transportation
%             - loading and cutting problems
%             - mining and process industries
%             - planning problems
%             - public services
%             - scheduling problems
%             - telecommunication
%             - time tabling
%
%         - A dynamic routine for generation of global test problems added (GLKS).
%
%        TOMLAB /SOL:
%
%         - SQOPT7 updated - an incorrect optimal value was reported for
%         some cases.
%         - Alpha release of SNOPT/SQOPT 7.1-1.
%         - QN (Quasi-Newton CG) QP Subsolver enabled.
%
%        TOMLAB /CGO:
%
%         - No longer estimating gradient at end of run.
%         - Random strategy updated.
%         - Updated behavior when replacing poor values.
%
%        TOMLAB /CONOPT:
%
%         - Automatically sets LS2PTJ depending on user information.
%         - Unused variables removed.
%         - Possible to set equality tolerance in Prob.CONOPT.eqTol.
%
%        TOMLAB /CPLEX:
%
%         - xState, bState reported correctly.
%
%        TOMLAB /MAD:
%
%         - isfinite, isnan, isvector and isinf added to distribution.
%         - Dummy files moved to mad/@dummy.
%         - tril, triu, cumsum added to supported functions.
%         - reshape, subsref updated.
%         - Support functions test_equal, zeroslike and isactive added.
%
%        TOMLAB /MINLP:
%
%         - Printout removed when optPar(20) was set for filterSQP.
%
%        TOMLAB /KNITRO:
%
%         - Help updated.
%
% 050922 maintenance release 4.9.1
%
%   TOMLAB /MAD:
%
%        * @fmad operations added: isfinite, isinf, isnan, isvector
%        * Dummy files are now in tomlab/mad/@dummy, to avoid problems if paths are set with subfolders added
%
%
% 050914 TOMLAB v4.9 released:
%
%        * glbDirect problem with NaN values fixed.
%        * Demo improvements for TOMLAB /LDO.
%        * One more MIQP problem added to distribution.
%        * One more LLS problem added to tomlab\testprob. LSQLIN fails on this case.
%        * cls_prob, con_prob updated.
%        * slsSolve much faster with SNOPT7 as a subsolver - CG used instead for optPar(66).
%        * Safeguarded starting point for PDCO modified.
%        * COLMMD replaced by COLAMD in ComputeQR, ssqls, cne, csne, colamd and qls.
%        * Example updated for expSolve.
%        * chs_prob updated (Hessian fixed for test case 104).
%        * isstr replaces by ischar in applicable files (30).
%        * Two new algorithms available in clsSolve - Li-Fukushima MBFGS and Broydens method.
%        * expSolve now takes eType as a new input.
%        * Improved error messages if problems with user routines.
%        * Weighting for least squares problems updated - numerical differentiation caused double weights.
%        * Automatic differentiation for least squares problems updated.
%        * Help updated for TLSQR.
%        * 'help solvername' will now display suitable information for most solvers.
%        * New quick guide released - detailed information about patterns, derivatives and warm start provided.
%        * BINTPROG added to optimization toolbox interface.
%        * L1LinSolve and L1Solve fully sparsified. Code more efficient now.
%        * infSolve - patterns now always sparse.
%        * goalSolve - improved error messages and sparsified code.
%        * slsSolve - code improved for speed.
%        * miqpAssign, miqqAssign and minlpAssign updated.
%        * Number of Trials in estConsPattern and estJacPattern reduced to 2.
%        * estd2cPattern added to distribution (needed when estimating d2LPattern).
%
%
%    TOMLAB /SOL:
%
%        * optPar 55-62 removed from the solvers. Restarts need to be done directly from TOMLAB.
%        * LSSOL printing problem removed for problems with no linear constraints.
%        * LPOPT/LP-MINOS/QP-MINOS/QPOPT/SNOPT/SNOPT7/SQOPT/SQOPT7 now accepting problems without linear constraints.
%        * StateDef revised for SNOPT 6.
%        * All version of SNOPT using nonderivative line search when Prob.CheckNaN is set.
%        * LSSOL, MINOS, LP-MINOS, QP-MINOS, SNOPT all have updated evaluation counters.
%
%
%    TOMLAB /NLPQL:
%
%        * Problem with demo license fixed.
%
%
%    TOMLAB /CPLEX /XA /CONOPT /Xpress:
%
%        * Try-catch statement added to identify installation problems.
%
%
%    TOMLAB /MINLP:
%
%        * Major update. Issue resolved regarding defaults (correct defaults provided at all times).
%        * filterSQP constraint evaluations corrected.
%        * Memory problem resolved for filterSQP (solving a small problem after a large could cause problems).
%
%
%    TOMLAB /PENOPT:
%
%        * Major update. Defaults now specified correctly.
%
%
%    TOMLAB /NPSOL:
%
%        * Return code revision. The solver sometimes returned feasibility for infeasible problems.
%
%
%    TOMLAB /MAD:
%
%        * isscalar added to functions supported.
%
%
%    TOMLAB /LGO:
%
%        * Default time limit changed to 1e7.
%
%
%    TOMLAB /CPLEX:
%
%        * Updated cpxtomtest1 to correct format.
%
%
%    TOMLAB /KNITRO:
%
%        * Updated handling of d2LPattern. Solver should perform better when all parts are analytical.
%        * d2LPattern automatically set for LP and QP problems.
%        * HessPattern automatically set for problems with linear or quadratic objectives.
%
%
%    TOMLAB /CONOPT:
%
%        * Hessian and d2c handling improved.
%        * Constraints only calculated if needed.
%        * Working correctly for all QP problems.
%        * No size restrictions to demo version anymore.
%        * d2LPattern automatically set for LP and QP problems.
%        * HessPattern automatically set for problems with linear or quadratic objectives.
%
%
% 050706 TOMLAB v4.8 released:
%
%         - SIGNIFICANT SPEED IMPROVEMENTS (TOMLAB /SOL) - solvers doing
%           (non-costly) callbacks to user routines now 20-50 % faster.
%         - preSolve updated to handle single constraints.
%         - M-file help improved system-wide.
%         - Manual more detailed about user supplied parameters.
%         - LPSOLVE renamed to MILPSOLVE.
%         - Optimization toolbox interface updated to handle more cases.
%         - Contents.m files now present in all folders.
%         - PrintResult working better when there are several solutions.
%         - DualSolve updated to handle changed version of lpSolve (now lpSimplex).
%         - Safe-guarded starting point added to mipSolve.
%         - Geometric programming added as a general problem type.
%         - As a result of Matlab Compiler issues, pragma %#mex added to all
%           solvers with m-file help for dll/mex.
%         - lpconAssign and qpconAssign modified to standard notation format.
%         - llsAssign safe-guarded for y and t input - they are now automatically
%           converted to full vectors.
%         - Geometric programming problem added to the TOMLAB Quickguide.
%         - MAD examples added to TOMLAB Quickguide.
%
%        TOMLAB /CGO:
%
%         - rbfSolve and ego has improved initial strategies.
%
%        TOMLAB /NPSOL:
%
%         - Specs file functionality updated.
%         - Minor bug related to LSSOL and NLSSOL fixed.
%
%        TOMLAB /SOL:
%
%         - No longer possible to use general TOMLAB parameters to set options.
%         - The optPar vector needs to be used for all alterations.
%         - Full manual review of control parameters.
%
%        TOMLAB /SNOPT:
%
%         - Return code and ExitFlags updated for SNOPT7 and SQOPT7.
%         - New parameters made available for SNOPT and SQOPT.
%         - Beta version of SNOPT 8 added to distribution.
%         - Several parameter and manual updates for all SNOPT version (6,7 and 8).
%         - SNOPT 6 has control vector (optPar) of length 65, SNOPT 7 of length 71.
%
%         - SNOPT 6 information:
%
%             Correction of defaults for optPar 2,10,27,28,29,34,35,36,42
%             New parameter LU SWAP TOLERANCE as optPar(25), LU DENSITY removed.
%             New parameter PROXIMAL POINT METHOD as optPar(64)
%             New parameter PENALTY PARAMETER as optPar(65)
%             New parameter NEW SUPERBASICS (MINOR SUPERBASICS) as optPar(66)
%             optPar(44): FEASIBLE EXIT (Obsolete)
%             Help added comments about parameters only available using SPECS file.
%
%         - SNOPT 7 information:
%
%             Correction of defaults for optPar 2,10,27,28,29,34,35,36,42
%             New 6.2 parameter LU SWAP TOLERANCE as optPar(25), skip LU DENSITY
%             New 6.2 parameter PROXIMAL POINT METHOD as optPar(64)
%             New 6.2 parameter PENALTY PARAMETER as optPar(65)
%             New 6.2 parameter NEW SUPERBASICS (MINOR SUPERBASICS) as optPar(66)
%             optPar(44): FEASIBLE EXIT (Obsolete)
%             Adding comments about parameters only available using SPECS
%             New 7.1 parameter: optPar(66). QPSOLVER CHOLESKY,CG or QNCG
%             New 7.1 parameters: optPar(67): CG TOLERANCE, optPar(68): CG ITERATIONS
%             New 7.1 parameter: HESSIAN (or CG) PRECONDITIONING, optPar(69).
%             New 7.1 parameter: SUBSPACE, optPar(70).
%             New 7.1 parameter: optPar(71): HESSIAN DIMENSION (also REDUCED HESSIAN)
%             Removed special setting of SUPERBASICS (optPar(48)).
%
%        TOMLAB /MINOS:
%
%         - Control vector (optPar) now of length 71.
%         - Added new parameters AIJ CONVERGENCE, SUBSPACE, LU WRAP.
%         - Change defaults for #6,42,43 (optPar control vector).
%         - Increased SUPERBASICS, optPar(48) to avoid too low value.
%         - LP-MINOS and QP-MINOS updated with new parameters.
%
%
%        TOMLAB /CPLEX:
%
%         - Version 9.1 released. See the manual for more information.
%         - Possible to warm start solution process for MILP and MIQP problems.
%         - New local branching heuristics.
%
%        TOMLAB /PENBMI /PENSDP:
%
%         - Warning messages removed.
%
%        TOMLAB /AMPL:
%
%         - Safe-guarded in Jacobian calculation.
%         - Both minimization and maximization problems possible to solve.
%
%        TOMLAB /KNITRO:
%
%         - Updated to choose appropriate HESSOPT depending on user routines
%           supplied.
%
%        TOMLAB /CPLEX:
%
%         - Issues warning of quadratic problem (hessian) is not symmetric.
%         - Possible to supply starting vector for mixed integer problems.
%         - Problems corrected for sparse matrix handling in MATLAB R12 and
%           R13 (only applies to problems with quadratic constraints).
%
%        TOMLAB /GP:
%
%         - Released
%         - GP is an interior-point package for geometric programming. Problems which are non-differentiable in the optimum can be efficiently solved.
%
%        TOMLAB /MINLP:
%
%         - Supported for Linux 64-bit and MAC OS X.
%
%        TOMLAB /CONOPT:
%
%         - Now supported for Sun Solaris.
%
%        TOMLAB /NLPQL:
%
%         - Now supported for all operating systems.
%
% 050601 Maintenance release - TOMLAB v4.7.2
%
%          - Minor bug in NPSOL mex file fixed.
%          - Final SNOPT 7 version released, minor changes from TOMLAB 4.7.0
%
% 050520 Maintenance release - TOMLAB v4.7.1
%
%          - A few missing files added.
%          - Corrupted MAT file fixed.
%
% 050509 TOMLAB v4.7 released:
%
%         - Example for calling TOMLAB solutions from Excel included.
%         - More standalone examples included in the distribution. An Excel
%           example is now included.
%         - Minor name changes to avoid conflicts with other toolboxes:
%              printmat -> PrintMatrix
%              assign -> ShowAssignment
%              simplex -> ShowSimplex
%              aic -> AkaikeIC
%              runfleq1 -> Trunfleq1
%         - New folder tomlab/common includes all files common to any installation.
%         - Only one startup file included in the distribution, all others
%           have solver specific names.
%         - Linear/Quadratic problem with nonlinear constraints added to quick guide.
%         - Simulation problem added to quick guide.
%         - Manual updated to clarify user parameters in Prob.
%         - uhs_prob, chs_prob and uc_prob updated.
%         - checkFuncs added. We recommend that this routine is used to check
%           TOMLAB problem before execution.
%              Prob = *Assign();
%              Prob.user.a = a;
%              ...
%              checkFuncs(Prob);
%              ...
%         - chs_prob, con_prob, glc_prob, minlp_prob linear constraints updated.
%         - Makeinitfile modified to avoid print-outs to screen.
%         - glcCluster updated to accept maxFunc3 input directly.
%         - glcFast updated to handle all variables fixed.
%         - New features for sending user parameters to Init Files, see Section
%           14.2 in the general TOMLAB manual.
%         - Compatibility verified with MATLAB R14 SP2.
%         - Review for MATLAB R12 - CPLEX, PENSDP, PENBMI working correctly.
%         - Problem with Prob.NumDiff = 4 resolved. The spline is now converted
%           correctly.
%         - A dedicated version for Linux on MATLAB R14 released.
%         - glcFast safe-guarded for user routines with different number of
%           inputs.
%         - lpSolve renamed to lpSimplex.
%         - LPSOLVE (LP and MILP solver) added to the TOMLAB Base Module.
%           See http://tomopt.com/products/base/solvers/LPSOLVE.php and
%           TOMLAB manual for more information.
%         - 19 more glb_prob's added to the distribution.
%         - New version of glcSolve released.
%         - glbDirect and glcDirect released (will replace glbFast and
%           glcFast eventually).
%         - Printing error in nlp_r corrected.
%         - L1LinSolve revised to avoid an intermediate LS structure
%           (improved memory handling).
%         - TOMLAB/MATLAB engine call from C program included in the distribution.
%         - Further devlopment of a constrained mixed-integer DIRECT
%           algorithm has been implemented in the new releases of glcSolve
%           and glcDirect.
%
%        TOMLAB /CGO:
%
%         - Updated to better handle return codes from glcFast.
%         - Improved handling of fixed variables.
%         - Possible to run 1-dimensional integer problems efficiently.
%         - Checks on duplicate points added.
%         - Stops if all integer combinations tried.
%
%        TOMLAB /SNOPT:
%
%         - Version 7 (7.1-1(5)) of SNOPT and SQOPT is now available.
%         - SNOPT 7 is called with tomRun('snopt7'...)
%         - Significant performance improvements - much larger problems can
%           now be solved.
%         - Documentation updated for scaling, 2 default for LP, 0 for NLP
%           (good to try scaling for NLP's manually)
%         - SQOPT documentation updated for scaling, 2 default for LP, 0 for
%           QP (good to try scaling for QP's manually)
%
%        TOMLAB /CPLEX:
%
%         - Examples aircrew, biptest, iptest, tomtest1, tomtest2 have new
%           names, cpx prepended.
%         - Documentation updated for sensitivity analysis - possible to set
%           indices directly rather than the start and end points.
%         - New options added to the solver, guided dives and more.
%
%        TOMLAB /Xpress:
%
%         - Examples aircrew, biptest, iptest, tomtest1, tomtest2 have new
%           names, xp prepended.
%
%        TOMLAB /AMPL:
%
%         - New manual available.
%
%        TOMLAB /PENSDP + /PENBMI:
%
%         - Now possible to turn off screen printing.
%
%        TOMLAB /OQNLP:
%
%         - Print out bug removed for feasibility tolerances.
%         - Works for pure integer/discrete programming problems.
%         - Feasibility problems solved.
%         - New manual available.
%         - Defaults now different for DISTANCE_FACTOR.
%         - More stable release - several issues fixed.
%
%        TOMLAB /MAD:
%
%         - Bug in multiplication removed.
%
%        TOMLAB /NPSOL:
%
%         - NLSSOL safe-guarded for incorrect length in residuals.
%         - Constraint handling updated for NLSSOL.
%
% 050202 TOMLAB Base Module v4.6 released:
%
%          - Bug in nlpSolve removed. Convergence testing was incorrect.
%          - Penalty values in conSolve now dampened to promote convergence.
%          - conSolve now works better for unconstrained problems.
%          - Major updates for goalSolve, several new subsolvers are now working.
%            Better handling when no analytical derivatives are available.
%          - glcCluster has slightly modified input fields and better defaults.
%          - Pattern checks added to checkDerivs.
%          - chs_prob set extensively corrected.
%          - Major revisions of L1Solve, infSolve, slsSolve and goalSolve.
%          - PDCO and PDSCO updated.
%          - WeightType 3 removed from LSEI.
%          - Tfzero and Tfmin conformed to TOMLAB standard.
%          - Complete print file review of the entire system. Minimal file
%            printing by default from all solvers.
%          - Variable, linear/nonlinear constraint states revised for all solvers.
%          - lpDemo and qpDemo updated to new TOMLAB format.
%          - Several updated to LDO to conform to TOMLAB structure.
%          - All assign routines does length checking on x_0, x_L and x_U.
%          - Starting point safe-guarded in all solver packages.
%          - tomSolve updated with more solver options.
%          - fmincon interface (drop-in replacement) updated to handle more cases.
%          - Numerical differentiation for NLLS problems improved.
%
%        TOMLAB /CPLEX:
%
%          - Dual CPU activated in demo license. Please contact
%            support@tomopt.com for more information.
%          - Possible to warm start the solution process for LP problems.
%
%        TOMLAB /PENBMI:
%
%          - Version 2.0 released.
%          - Safeguarded for missing initial vector.
%
%        TOMLAB /PENSDP:
%
%          - Version 2.0 released.
%
%        TOMLAB /SOL:
%
%          - New manuals for most solvers (MINOS, SNOPT, QPOPT, SQOPT).
%
%        TOMLAB /SNOPT:
%
%          - Default for number of superbasics is automatically set to
%            better number.
%          - Safeguard added in SQOPT for dense QP problems.
%
%        TOMLAB /NPSOL:
%
%          - WeightType 3 removed from LSSOL.
%          - Print files secured in the code to avoid unexpected behavior.
%
%        TOMLAB /KNITRO:
%
%          - Minor interface upgrades.
%          - Stability improvements.
%
%        TOMLAB /NLPQL:
%
%          - Print files no longer generated by default.
%
%        TOMLAB /OQNLP:
%
%          - Version 3.0 released.
%          - Direct access to MSNLP from the interface.
%          - The embedded local solver LSGRG-2 can be executed independently.
%          - Empty print file removed.
%          - More detailed algorithmic description available in the TOMLAB /OQNLP manual.
%          - Possible to generate a LOCALS file with all local solutions.
%          - Return codes reviewed.
%
%        TOMLAB /MINLP:
%
%          - WarmDefDUNDEE updated with estimated problem patterns.
%          - New definition of infinity for better performance.
%          - New manual with algorithmic description available for download.
%
%        TOMLAB /XA:
%
%          - Variable and linear constraint states added to output.
%
%        TOMLAB /LGO:
%
%          - New version with improved stability.
%          - Solver options are now set in Prob.LGO.options.
%
%        TOMLAB /CONOPT:
%
%          - Solver options are now set in Prob.CONOPT.options.
%
% 041130 TOMLAB v4.5 released. News:
%
%        TOMLAB Base Module:
%          - WarmDefDUNDEE added for the TOMLAB /MINLP solvers.
%          - Direct access to the tomsol.dll enabled - improving speed for
%            matrix - vector mutliplications in MATLAB/TOMLAB. tomsol(0) will
%            display the help.
%          - PrintAssign working with user given names for problems.
%          - 30 additional test cases for MILP problems included in the
%            distribution. See mip_prob.m
%          - 2 additional MIQP problems added.
%          - 42 additional LP problems added.
%          - 27 extra QP problems added.
%          - lpconAssign and qpconAssign added for LP and QP problems with
%          - nonlinear constraints.
%          - xnargin modified provide better error messages.
%          - xxx_prob.m removed from distribution and integrated in applicable files.
%          - Bugs in estHessPattern/estConsPattern/estJacPattern fixed.
%          - fmincon interface updated.
%          - tomRun changed to tomRun('solver', Prob, PriLev), i.e. no
%            need to give a third [] parameter.
%          - TLSQR updated to exit more smoothly on linear least squares
%            problems with equality constraints.
%          - mipSolve and cutplane can now solve LP problems.
%          - WARNING - glcAssign has been changed, setupFile, nProblem and
%            KNAPSACK has been removed as inputs - WARNING
%
%        MATLAB Compiler:
%          - C/C++ examples for Windows/Linux/Sun included in the general
%            distribution. This examplifies how to call TOMLAB in standalone
%            mode.
%          - An updated TOMLAB /SAL guide is available from the download page.
%
%        TOMLAB /KNITRO v4.0 released:
%          - Many new features, including a new solver algorithm option.
%          - Defaults changes to better suit the user. If analytical gradient,
%            then finite-difference Hessian-vector products are used. If only
%            the objective function is given, a memory limited quasi-Newton
%            Hessian is used.
%          - A few minor fixes included in the new version.
%          - Download the latest manual for more information.
%
%        TOMLAB /SOL:
%          - New defaults for line search when using numerical derivatives.
%          - Updates to handling of defaults for the solver settings.
%
%        TOMLAB /CGO - Version 2.5 released:
%          - rbfSolve updated for better handling of binary and integer variables.
%          - EGO improved to work better with bad input data.
%          - New manual with detailed algorihtmic descriptions for rbfSolve.
%          - A problem with rbfSolve and ego sampling outside the bounds for
%            costly integer programming has been removed.
%          - Safeguards added when local solver cannot return integer solution.
%
%        TOMLAB /XA v14 released:
%          - New features for IIS.
%          - Barrier works for LP problems.
%          - Faster simplex and MIP solutions.
%
%        TOMLAB /LGO
%          - MaxCPU added. Minor printing bugs removed.
%
%        TOMLAB /CPLEX and TOMLAB /Xpress
%          - (X)MPS and other input files can now be read into MATLAB.
%            LP, MILP, QP and MIQP supported.
%          - The AMPL format is supported in TOMLAB /AMPL.
%
%        TOMLAB /MAD
%          - New version of mtimes_dab. ceil and cumprod added to the distribution.
%
%        TOMLAB /AMPL
%          - Printlevel added to amplAssign.
%          - Safeguard for cases with no nonlinear constraints.
%
%        TOMLAB /OQNLP
%          - New interface routine included which avoids numerical
%            differentiation for integer variables.
%
%        TOMLAB /MINLP
%          - Minor printing fixes.
%
%        TOMLAB /CONOPT
%          - Released for Linux.
%
%        TOMLAB no longer supported for MAC OS 9. Existing licenses can be
%        transfered to MAC OS X at no cost.
%
% 040929 Maintenance release. Changes:
%
%        MINOS, NPSOL, NLSSOL, LSSOL restore default solver parameters on
%        consecutive runs.
%
% 040928 TOMLAB v4.4 released. News:
%
%        Infeasibility and sensitivity analysis in TOMLAB /CPLEX
%        and TOMLAB /Xpress.
%
%        Tomlab /MAD v1.2 - general improvements, recodes of certain
%        functions
%
%        Matlab Compiler 4 compatibility improved.
%
%        Minor bug fixes to TOMLAB /CGO, /XA, /XPRESS, /CPLEX
%
% 040828 TOMLAB CPLEX 9.0.2. Incompatibilites with Matlab R14 fixed.
%
%        TOMLAB /CPLEX new features: Sensitivity and Infeasibility Analysis
%        Minor problems with TOMLAB under Matlab 7 have been fixed.
%
% 040728 TOMLAB SAL modification
%
%        Several pragmas and file changes in TOMLAB to better support the
%        MATLAB Compiler v4. mFiles now tomFiles.
%
% 040521 TOMLAB /OQNLP v2.0 released.
%
%        Many general improvements to the multi-start features. Changes
%        to the distance and merit filters make them theoretically sounder
%        and more adaptive to the problem instance. Stochastic drivers are
%        implemented whose performance is comparable to the OptQuest scatter
%        search implementation, if reasonably tight bounds on all variables
%        are imposed.
%
% 040428 TOMLAB /PATH released.
%
%        Linear and nonlinear mixed complementarity problems can now be
%        solved with TOMLAB. The new assign routines lcpAssign and mcpAssign
%        should be used to create the problem.
%        TOMLAB /PATH also handles linear and convex quadratic programming
%        problems. Developed in cooperation with University of Wisconsin
%        at Madison.
%
% 040417 TOMLAB v4.3 released.
%
%        The TOMLAB Base Module now features several new additions for
%        derivative handling:
%        Automatic estimation of ConsPattern, HessPattern and JacPattern for
%        large-scale problems.
%        New method for numerical differentiation using standard MATLAB splines.
%        The TOMLAB interface routines now supply the user with row, column and
%        variable information for use in their code. With this information, only
%        segments of costly code need to be computed.
%        Many new feature does automatic problem validation. Element sizes in
%        linear constraints can be checked, the PreSolve capabilities have been
%        improved, a special flag for warning messages is now used.
%        A MaxCPU flag has been introduced for several solvers, limiting the
%        amount of CPU time consumed before returning to the command line.
%
%        TOMLAB /CPLEX has been updated with a special network interface.
%        The user can now specify a set of nodes and arcs when solving their problems.
%
% 040325 TOMLAB released for MAC OS X.
%
% 040323 TOMLAB /CGO v2.0 released.
%
%        The CGO solvers, rbfSolve and EGO now supports integer variables.
%        TOMLAB /OQNLP can be used as a sub solver.
%
%        TOMLAB folder structure changed, new directories: base and cgo
%
% 040307 TOMLAB /MAD released.
%
%        MAD is not available for all users in standalone mode. The TOMLAB
%        Base Module features a complete integration with MAD by the use of
%        flags in the Prob structure.
%
% 040226 TOMLAB /NLPQL released.
%
%        NLPQLP, NLPJOB and DFNLP included. NLPQLP is an SQP solver
%        NLPJOB is a package for multi-criteria optimization, while DFNLP
%        is mainly for nonlinear fitting problems.
%
% 040204 Release of Tomlab 4.2
%
%        Now possible to call MEX solvers directly through TL files,
%        provided the Prob structure has been ProbCheck'ed. Example:
%
%           Prob   = ProbCheck(Prob,'snopt');
%           Result = snoptTL(Prob);
%
%        SNOPT and SQOPT version 6.2-2 (previously 6.1-1)
%
%        LSSOL minor fix in MEX file, could potentially cause crashes
%
%        Tlsqr bug in PDCO mode fixed.
%
%        Tomlab now ready for integration with MAD Toolbox for Automatic
%        Differentiation
%
%        Bug fixes in TOMLAB /XA.
%
%        Upgraded versions of glbSolve, glbFast, glcSolve, glcFast,
%        glcCluster.
%
%        TOMLAB /MINOS now support LU Rook Pivoting.
%        TOMLAB /SNOPT now support LU Rook Pivoting, and LU Diagonal Pivoting.
%
%        TOMLAB /CGO v1.6 released. Better handling of non-costly nonlinear
%        constraints implemented.
%
%        New version of TOMLAB /CPLEX v9.0. Explicit quadratic constraints
%        are now supported. With a special license the solver may be excuted
%        on up to 64 parallel processors.
%
%
% 040128 Release of TOMLAB /LGO.
%        Global optimization package.
%
% 031201 Release of 4.1.2
%        Inf on variables could cause problems for miqpBB
%
% 031027 Release of 4.1.1
%
%        tomHelp system rewritten.
%
%        tomRemote feature for running with a menu-style over text-only connections.
%
%        clsSolve: spdiags used for QP subproblems when >1000 variables.
%
%        FDJac.m : sends index of perturbed variable for possible speedup in
%                  calculating numerical constraint derivatives.
%
% 030905 Release of 4.1.0
%
%        New solvers added: CONOPT, KNITRO, OQNLP
%
%        mipSolve: Dual gap now handled correctly.
%
% 030526 Release of 4.0.6
%
% 030525 simAssign,sim_fc,sim_gdc,sim_f,sim_g,sim_c,sim_dc: Special handling
%        of e.g. simulation problems, where both f and c must be computed
%        at the same time.
%        Change ProbDef, ProbCheck, mFiles,
%        Change many xxxAssign routines to use mFiles to set user files.
%
% 030525 clsSolve: Add two more algorithms, not fully tested
%
% 030525 glcPrint, glcFast: Add print of all sampled x, revise print levels
%
% 030524 rbfSolve, rbf_xxx: Revision of constraint handling
%
% 030514 ego: Revision of constraint handling
%
% 030427 Add solver NLPQL
%
% 030228 Release of 4.0.5
%
% 030227 TPvogel. Failed computing bfs for TPsimplx. 2 lines were commented
%
% 030226 rbfSolve: Minor errors corrected, change of variables names.
%
% 030221 Release of 4.0.4
%
% 030221 Tlsqr: nOut undefined
% 030221 rbfSolve: . missing
%
% 030220 minlpBBTL: New input parameters mlp and maxf included
% 030220 qp_Hess: Illegal character in file
% 030220 nlp_d2c: Avoid unnecessary computations
% 030220 qp_Hess: Illegal character in file
%
% 030214 Release of 4.0.3
%
% 030211 Adding Hoch-Schittkowski unconstrained and constrained test problems
%        In uhs_prob and chs_prob
%
% 030211 nlp2_fgH, ls2_rJ, ls2_rJS: Updated to v4.0, using global NARG
%
% 030211 PrintResult. Avoid printing Hessian information, when not used
%
% 030211 lseiTL. Avoid printing Hessian information, when not used
%
% 030210 Release of 4.0.2
%
% 030210 Update of penfeas_lmi and penfeas_bmi in /PENSDP, /PENBMI
%
% 030208 Two new MINLP test problems, minlpDemo, minlpAssign added.
%
% 030206 /MINLP: Use filterSQP as solver name, correct comments in all /MINLP
%        related files. Change names to filterSQPTL.m, filterSQP.m
%
% 030206 FDHess: Zero matrices not correctly treated.
%
% 030203 Release of 4.0.1. Minor bug fixes.
%
% 030203 TomlabVersion: Vector too short for old licenses. Always expand with 0s
%
% 030131 Make Tomlab handle function_handle function input. Changes only
%        in xnargin.m and xxnargin.m.
%
% 030131 fdng3: Loop variabel i in conflict with complex i
%        fdng: In some cases g got to be a vector. Initialize as column of 0s.
%
% 030131 Release of 4.0.0, /MINLP v1.0, PENBMI v1.0, PENSDP v1.0.1
%
% 030129 Make Optimization Toolbox interfaces compatible. Handle matrix x,
%        function handles, inline expression, @ functions, cell arrays.
%
% 030127 Use global structure otxProb to send Tomlab structure information
%        to Optimization Toolbox solver interfaces. Use Prob.Solver.Tomlab
%        to change the default solver used. Still possible to send structure
%        as first extra argument to solver as well.
%
% 030127 If Prob.CheckNaN ~=0, nlp_d2c, nlp_H, nlp_d2r checks for NaN elements
%        and estimates the corresponding derivatives numerically.
%        If Prob.CheckNaN >0, the same applies for nlp_dc, nlp_g, nlp_J
%        Off-diagonal elements in symmetric Hessians should both be set as NaN.
%        New versions of fdng, fdng2, fdng3, estimate only NaN elements in
%        gradient, if gradient vector is input.
%
% 030118 Change name of dfzero to Tfzero
%
% 030117 Major revision for v4.0. Integrate /PENBMI
%
% 010117 If CheckNaN==1, check for NaN in derivatives, and estimate any
%        NaN elements numerically
%
% 030117 CreateTomProb: Use probtype 14 as Linear SDP with BMI, not MIQQ
%
% 030116 Change names: fmin to Tfmin, wnnls to Tnnls, lsqr to Tlsqr
%
% 030113 nlp_d2c, FDcHess: Estimate constraint Hessian numerically
%
% 030110 tomRun: Call preSolve if Prob.optParam.PreSolve set
%
% 030107 mipSolve: Avoid large arrays being allocated
% 030107 mipSolve: Empty b_L gives b_L=b_U, correct comments
%
% 021217 L1LinSolve ...
%
% 021021 Release of CGO 1.4
%
% 021020 ego: Having linear, but not nonlinear constraints caused crash, fixed.
%        Revision to handle infeasibility. Include Gutmann initial strategy.
%        ego_cc: Must test on missing nonlinear constraints.
%
% 021020 rbfSolve: Not setting Prob.GO caused crash, fixed.
%        Revision to handle infeasibility
%
% 021011 Release of 3.2.2 and CGO 1.3
%
% 021010 nlpgui3: Avoid help menus on MAC, not working in Matlab 5.2
%
% 021010 checkType: Change of logic, otherwise causing crash in tomGUI.
%
% 021001 optim_r: Reorder of statements, otherwise information overwritten
%
% 020929 fmincon: Add constraint tolerance TolCon - Prob.optParam.cTol
%
% 020928 goalSolve,goal_c, goal_dc: New multi-objective goal attainment routine
%
% 020921 Preparation for CPLEX in Tomlab
%
% 020920 defblbu: Using wrong dimension in some cases
%
% 020919 nlp_cdcS: Check on Prob.ConsDiff == 6, avoid Tomlab differences
%
% 020901 Release of CGO 1.2.
%
% 020831 rbfSolve: Major revision. Input and search rule extensions
%
% 020831 SolverList: Major revision. Gives list of licensed and nonlicensed
%        solvers for each problem type
%
% 020830 Release of 3.2.1
%
% 020828 rbfSolve, ego: Use GetSolver to select default local solver
%
% 020824 rbfSolve, ego: Only accept local solutions if linear constraints OK
%
% 020823 rbfSolve, ego: Error transforming linear constraints if SCALE=1
%
% 020821 quadprog. Financial TB sends empty tolerances (options.TolX).
%        Check first if empty.
%
% 020821 rbfSolve, ego: Local surface search from all points (max 20)  found
%        by global search, not just first point. Change of amount of printing
%
% 020820 lpSolve, conSolve, qpSolve, nlpSolve, ucSolve, clsSolve, systest,
%        StateDef: Changes due to Matlab 6.5 logical handling
%
% 020705 Release of 3.2.0 and CGO 1.1. New license handling
%
% 020701 conSolve. Use full(B_k) in call to pinv
%
% 020701 GetSolver revised, more general solver selections
%
% 020701 Add problem types miqp, miqpp, minlp & sdp. Changes in GUI, mat-files
%
% 020701 checkdll, tomlabVersion revised for new license handling
%
% 020630 For QPs with large, dense or nearly filled quadratic F matrices,
%        it is faster to make a callback from the MEX to compute F*x in Matlab.
%        The F matrix is then never copied to the MEX.  Prob.SOL.callback = 1
%        for SQOPT, and Prob.DUNDEE.callback = 1 enables this feature.
%
% 020630 Tomlab now prepared for solvers filterSQP, minlpqq, miqpbb and PENSDP
%
% 020628 nlp_f,nlp_g, nlp_H, iniSolve.m. More efficient handling of separable
%        functions. Fixed bugs creating unnecessary calls. Speedup.
%
% 020616 nlp_fg. Avoid call to nlp_g if NumDiff = 6, check was removed in nlp_g
%
% 020616 fdng. Test to avoid division with 0 for sick calls.
%
% 020610 rbfSolve. Correcting errors in print statements. Avoid square root of
%        negative numbers when having ill-conditioning (in tomsol MEX)
%
% 020613 iniSolve, nlp_f, nlp_g, nlp_H. Use global variable PartSep to speed
%        up handling of separable problems
%
% 020531 clsSolve. Set LargeScale=0 if n<10 to avoid bug in sqr2 package
%
% 020512 ProbCheck,ProbDef: Adding field Prob.DUNDEE, Fletcher/Leyffer solvers
%
% 020512 tomRun, mexRun, optParamDef: Adding entry for Dundee QP solver bqpd
%
% 020506 Release of 3.1.3
%        New versions of SNOPT, SQOPT, TOMSOL
%
% 020423 plotUtil. Updates to plot rbfSolve and ego sampled points. Special
%        plot for the initial points.
%
% 020422 glcSolve, glcFast. Better way to compute rectangle size for integers.
%
% 020422 Add expDemo. Exponential fitting with TQ and IF format.
%        Improved expLinCon, exp_prob and expInit.
%
% 020417 Add sqr2 sparse QR package. First use in clsSolve, that calls new
%        routine ssqls that is using sqr2. Avoids forming of huge orthogonal
%        Q matrix of size mxm, where m is number of residuals.
%
% 020417 PrintResult: Print index of worst constraint violation
%
% 020416 Revision in handling of numerical derivatives for fmincon interface.
%        Making new test example in testfmincon, new file fmincon_c.
%
% 020416 Revision in handling of numerical derivatives for infSolve, L1Solve
%        and slsSolve. Changes in FD*JAC.m, fd*ng.m, *=1,2,3
%        Change ?_dc, for ?=mmx,L1,L2, Change ?Solve, ?Demo for ?=inf,L1,sls
%
% 020414 Release of 3.1.2
%
% 020413 slsSolve, L2_???.m: New Sparse Least Squares solver.
%
% 020413 L1Solve, L1_???.m: New L1 solver.
%
% 020413 infSolve, mmx_???.m: Revised to handle nlp_c,nlp_dc calls correctly
%
% 020412 PrintResult: Use multiple lines for multiple global solutions
%
% 020411 glcCluster: Remove memory and speed bottlenecks with Fortran code
%
% 020409 tomlabInit, globalSave, globalGet, iniSolve. Use global variable
%        NARG to save xnargin("user function"), speeds up low level routines
%
% 020409 Speed up all low level routines nlp_f, nlp_c etc.
%
% 020409 LineSearch, xxxSolve, xxx=con,uc,cls,nlp: Set Prob.Mode, Prob.nState
%
% 020409 nlp_cdcS,nlp_cdc,nlp_fg: Send Prob.Mode, Prob.nState to low level
%        user routines, generalize Tomlab to handle mode similar to SOL solvers
%
% 020404 ego_cc, new routine to handle constraints as well as variable scaling
%        rbfSolve, ego: Revision for constraint handling
%
% 020404 infSolve, mmx_?: Major revision of infSolve for minimax problems.
%
% 020304 CreateProbQP: Set max iterations to max of user given, and default.
%        Sometimes default was too low when running many LPs in mipSolve.
%
% 020304 mipSolve: Avoid using LPOPT for LP sub problems, always use MINOS.
%        Wrong field name used made SOL.optPar not correctly set
%        Print Inform from solver.
%
% 020304 snoptTL, minosTL, npsolTL, nlssolTL: Set DERLVL dependent on
%        Prob.ConsDiff and Prob.NumDiff.
%
% 020228 sqopt: Generalized to handle implicitly defined QP problems, where the
%        quadratic term x'*H*x is returned as z = H*x in a callback routine.
%        Function HxFunc.m is an example callback routine, and the call to
%        sqopt is in comments in the sqoptTL.m routine. The only change is that
%        'HxFunc' is 4th input argument instead of Prob.QP.F.
%
% 020210 nlp_cdcS: More efficient and safe way to convert dynamic sparse
%        Matlab format to static Fortran format. Using precomputed Prob.ConsIdx
%
% 020208 WarmDefSOL. Spelling error, cLambda should be cLamda
%
% 020210 minosTL,snoptTL: Compute linear index for nlp_cdcS in Prob.ConsIdx
%
% 020204 mipSolve: Pick up Prob.SOL.optPar/PrintFile and use if SOL solvers
%
% 020113 Release of 3.1.1.
%
% 020113 Syntax error in DefineGUIMenus.
%
% 020113 New test routine CGOrun in examples.
%
% 020111 Release of 3.1.0, and new toolbox Tomlab /CGO 1.0
%
% 020111 PrintResult: Use exponential format for f(x) if |f(x)| > 1E30
%
% 020111 glcCluster: Safe guard to odd number input, print improved f(x) info.
%
% 020110 inisolve: Avoid building search directions if not GUI or MENU system
%
% 020110 nlp_f: Set global PBUILD=0 to avoid wrong calls from nlp_r to pbuild
%
% 020110 glbFast, glcFast: Change fMin to glbfMin, glcfMin. Conflict with
%        Matlab function fmin.m when running glcCluster and reading fMin.
%        glbSolve: Change fMin to glbfMin; glcSolve: Change fMin to glcfMin
%        glcCluster: Change fMin to glcfMin, Loop to get feasible in 1st step
%
% 020104 ego_c, dace_f: Speedup for new toolbox Tomlab /CGO
%
% 020104 fminsearch: Prob was wrongly picked from input.
%
% 020103 ego: Major revision and speedup. Made similar to rbfSolve
%
% 011228 tomGUI.m: Check output from findstr with isempty for AD check
%
% 011226 glcAssign, mipAssign, probCheck, mip_prob, abc2gap: Change fields
%        SOS1 and SOS2 to sos1 and sos2, to conform with Tomlab / Xpress
%
% 011214 Change to use tomlablic instead of tomlab as license file
%        Made a new tomlab.m that calls tomlablic
%
% 011213 mipAssign, glcAssign, ProbCheck. Field MIP.SC changed into two fields:
%        SC semi-continuous variables,0 and real interval, either 0 or [a,b])
%        SI semi-integer variables with integer interval, 0 or a,a+1,...,b
%        These type of variables are only handled by Tomlab /Xpress v2.0
%
% 011212 Avoid use of local variable ProbName in SOL interfaces.
%        glcSolve.m: Wrong variable name in never executed error part.
%
% 011212 Remove script files tomlab\lib\inibuild.m and tomlab\lib\endbuild.m,
%        they are now included in iniSolve.m and endSolve.m.
%
% 011206 Release 3.0.25
%
% 011206 conDemo: Added examples for NPSOL and SNOPT.
%
% 011206 rbfSolve: Update of cycle strategy, minor revision.
%
% 011206 conAssign: Default 0 for fLowBnd for line search, changed to -realmax.
%        Caused failure for conSolve on problems with f(x) < 0.
%
% 011205 exp_q.m: Wrong matrix size estimating initial values for 4 exp-terms
%        Safeguarding for eType == 4. Revision of exponential fitting in
%        exp_prob, exp_r, exp_J, could crash for some cases.
%
% 011205 LineSearch.m: Did not return all info when reaching max iterations
%
% 011204 nlp_fc, nlp_gdc. Change interface because of bug/limitation in constr
%        nlp_cF had bug handling constraint gradient interface to constr
%
% 011204 Revision of Opt tbx 2.x interfaces lsqlin,lsqnonlin
%        Added lsqlin to GUI and interface opt20Run.
%
% 011203 fmincon,fminunc,quadprog,linprog,fminsearch:
%        Revision of Optimization toolbox 2.x interfaces. Name conflict
%        in Prob structure. fmincon and fminunc interfaces changed to use field
%        Prob.Solver.TOM for the Tomlab solver name.
%        opt20Run.m: Use Tomlab Result output to get nice GUI printing of result
%
% 011130 tomGUI: Help did not work for all problems, corrected.
%
% 011130 mipProbSet: Wrong computation of Prob.N, number of variables.
%        Also safeguarded some of the other similar routines.
%
% 011114 Release 3.0.24
%
% 011112 glcCluster.m: New hybrid global optimization algorithm
%        rbfSolve.m: 1st version of global solver for costly f(x).
%
% 011112 ProbCheck.m: A rare bug fixed, setting Prob.N = 0 (length(x))
%
% 011109 Fixed bug in all MEX solvers. If a tomlab directory was in the current
%        directory, the tomlab license was not found
%
% 011107 New solver balas01.m. Solves 0/1-problems with Balas algorithm.
%        Calls ldo\balas.m (totally revised).
%
% 011031 Added maxTri (maximum rectangle size) in Result structure, and in
%        printing in PrintResult and the GUI.
%        nlpsub.m: Add input Solver to plotUtil for global optimization plots
%        plotUtil.m: Using input Solver determines mat-file w global opt-info.
%        Major revision of glbSolve, glcFast, glcSolve, glcFast.
%
% 011030 Made many test problem names shorter, making GUI readable on Linux
%        Made GUI wider, to handle longer strings in menus, especially Linux
%
% 010903 Release 3.0.23
%
% 010903 Now 63 parameters in Prob.SOL.optPar to allow for the new option
%        LU Complete Pivoting in snopt, minos, sqopt.
%        Thus Prob.SOL.optParN = 63 in CreateProbQP.m, ProbCheck.m, ProbDef.m
%        Additional default parameter defined in SOLGet.m
%        tomGUI: GUI changed to use LU Complete Pivoting parameter
%
% 010903 New versions of snopt, minos and sqopt. Amount of memory sometimes
%        too low for these solvers. Now fixed. By default all three are
%        now silent and do not create any log files.
%        snopt and minos now use size(Prob.ConsPattern,2) to determine number
%        of nonlinear variables in constrained Jacobian. This size may differ
%        from the total number of variables, i.e. no need to send large 0-block
%        Cold start changes of hs and nS OK for minos (NLP,LP,QP), snopt, sqopt
%
% 010903 DefineGUImenus.m: Added glcSolve,glcFast in GUI for optType 9, glb.
%
% 010815 glcFast.m: New solver glcFast added, fast MEX version of glcSolve.
%
% 010728 Release 3.0.22
%
% 010726 glbSolve.m: Minor speedups.
%
% 010726 startup.m: Now tomlab main directory first in path. Otherwise e.g.
%        "help contents" would show contents.m in admat directory.
%        An entry "if 0" for TOMLAB / Xpress added. Set as "if 1" if
%        startup should generate a path to the TOMLAB / Xpress toolbox.
%        Must be edited if path other than "D:\tomlab\xpress".
%
% 010726 GetSolver.m: Differentiate between sparse and dense NLP (con) problems
%
% 010726 ego.m: Change local / global solver handling and printing.
%
% 010725 glb_prob.m: Correction of many fGoal for problems 13-24.
%
% 010725 nlp_cdcS.m: Sparse matrix handling corrected when user is using the
%        the non-Tomlab format.
%
% 010725 snopt.m: Changed instructions for non-Tomlab format.
%
% 010718 ego.m: Changed to use a general global optimization solver set as
%        globSolve (e.g. glbSolve or glbFast). Added a call to a local solver
%        (npsol) to improve the global optimum found.
%        Use convergence test, checking fGoal, similar to glbSolve and glbFast.
%        Include code for dace_bounds.m as local function in ego.m.
%        Changed call to dace_init.m to new function daceInit.m, which gives
%        Latin square initial set of points. The function value for these
%        points now computed in ego after the call, to conform with rbfSolve.
%        Added comments about input / output parameters
%
% 010718 dace_init.m: Removed, now daceInit.m is used.
%
% 010718 daceInit.m: Made from dace_init.m. Only generates initial points,
%        does not compute the corresponding function values.
%        Lower and upper bounds now input, not the full Prob structure.
%
% 010718 ego_c.m: Add test if condition value is inf, then set a big number.
%
% 010717 tomSolve.m: Added solvers glbFast and xpress-mp.
%
% 010717 tomGUI.m, tomRun.m: Added solver glbFast.
%
% 010715 glbFast.m: New solver glbFast added, fast MEX version of glbSolve.
%
% 010715 glbSolve.m: Variable name changes to conform with Fortran names in
%        glbFast. The new names then used in glbSave.mat, used for warm starts
%
% 010715 optParamDef.m, DefineGUIMenus.m, SolverList.m: Added solver glbFast.
%
% 010528 qpAssign.m, lpAssign.m, mipAssign.m: The name of the problem was not
%        defined in Prob.Name when no Init File was created, now added.
%
% 010621 Release 3.0.21
%
% 010419 MAJOR RELEASE OF TOMLAB /SOL v3.0
