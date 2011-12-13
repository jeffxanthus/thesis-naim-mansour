% tomlabInit.m 
%
% Call: 
%
%     tomlabInit 
%
% before running the TOMLAB routines
%
% -------------------------------------------------------------------------
% Change the following variables to suit your own configuration:
%
% MAXCOLS  Number of screen columns. Default 120.
% MAXMENU  Number of menu items showed on one screen. Default 50.
% -------------------------------------------------------------------------
%
% Global counters: 
% n_f      The number of function evaluations
% n_g      The number of gradient evaluations
% n_H      The number of Hessian evaluations
% n_c      The number of constraint evaluations
% n_dc     The number of constraint gradient evaluations
% n_d2c    The number of constraint Hessian evaluations
% n_r      The number of residual evaluations
% n_J      The number of Jacobian evaluations
% n_d2r    The number of evals of 2nd derivative of residual

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Aug 17, 1997.   Last modified June 28, 2002.

% Here we declare the global variables used by TOMLAB 
% -----------------------------------------------------------------
global MAXCOLS MAXMENU
MAXCOLS=120;   % SET MAXCOLS to 120 if you have a large screen
MAXMENU=50;    % Used by strmenu if many items, like for CUTE test problems

% -----------------------------------------------------------------

global glbSave   % Used to save globals in recursive calls to TOMLAB
                 % See routines globalSave and globalGet
global GlobalLevel

% -----------------------------------------------------------------
% Counters of call to m-files
global n_f n_g n_H n_c n_dc n_d2c n_r n_J n_d2r

% -----------------------------------------------------------------
% Save result of xnargin in NARG - speed up low level routines
global NARG
% -----------------------------------------------------------------
% User communication between f,g,H,c and dc routines 
global US_A US_B 
% Communication nlp_f/g and nlp_c/dc
global NLP_x NLP_f NLP_xc NLP_c NLP_pSepIndex
% Used by ls_f, ls_g, ls_H, r/J comp routines 
global LS_x LS_r LS_xJ LS_J LS_A wLS
global SEP_z SEP_Jz
global NLP_xg NLP_g NLP_xH NLP_H

% -----------------------------------------------------------------
global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print
MAX_x=20;
MAX_c=20;
MAX_r=30;

% -----------------------------------------------------------------
% Used by pbuild to build search directions and line search steps from x points
global BUILDP p_dx alphaV X_min X_max F_X X_OLD X_NEW pLen

% -----------------------------------------------------------------
% optType  defines type of optimization problem, 
% probType defines the problem type, same as optType, or a simpler problem
% solvType defines solver type.
global optType probType
optType = []; probType=[];  % Always reset to non defined

% Used to store results when running GUI
global ResultGUI  

global runNumber plotData           % Used in GUI, tomGUI.m
global question answer instruction  % Internal use in GUI, tomGUI.m

% -----------------------------------------------------------------
% Other globals, used internally between TOMLAB routines
% Robot problem
% global Baux B1 B2 B3 B4 B5 B6 mass_vector_aux gravity_aux mass_matrix gravity

% -----------------------------------------------------------------
global NTS_Zt NTS_alpha NTS_yW  % Used in NTS