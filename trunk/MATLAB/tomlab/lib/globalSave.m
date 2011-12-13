% globalSave.m
%
% function globalSave(depth)
%
% Saves all global variables in a structure glbSave(depth) and then
% initialize all of of them as empty. To be used in recursive optimization,
% calling several instances of TOMLAB.
%
% By using a depth variable, an arbitrarily number of recursions are
% possible.
%
% globalGet(depth) retrieves the variables

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2006 by Tomlab Optimization Inc., $Release: 5.7.0$
% Written Nov 8, 1998.       Last modified Dec 12, 2006.

function globalSave(depth)

if nargin < 1
   depth=1;
end

global glbSave 

d=depth;

% -----------------------------------------------------------------
% Result of xnargin in NARG, of nargout in NARGO
global NARG NARGO
glbSave(d).NARG         = NARG;
glbSave(d).NARGO        = NARGO;

% -----------------------------------------------------------------
% User communication f,g,H,c,dc
global US_A US_B
glbSave(d).US_A         = US_A;
glbSave(d).US_B         = US_B;

% -----------------------------------------------------------------
% Communication nlp_f/g and nlp_c/dc
global NLP_x NLP_f NLP_xc NLP_c NLP_pSepIndex

glbSave(d).NLP_x        = NLP_x;
glbSave(d).NLP_f        = NLP_f;
glbSave(d).NLP_xc       = NLP_xc;
glbSave(d).NLP_c        = NLP_c;
glbSave(d).NLP_pSepIndex= NLP_pSepIndex;

NLP_x=[]; NLP_f=[]; NLP_xc=[]; NLP_c=[]; NLP_pSepIndex=[];

% -----------------------------------------------------------------
% Used by ls_f, ls_g, ls_H, r/J comp routines 
global LS_x LS_r LS_xJ LS_J LS_A wLS

glbSave(d).LS_x        = LS_x;
glbSave(d).LS_r        = LS_r;
glbSave(d).LS_xJ       = LS_xJ;
glbSave(d).LS_J        = LS_J;
glbSave(d).LS_A        = LS_A;
glbSave(d).wLS         = wLS;

LS_x=[]; LS_r=[]; LS_xJ=[]; LS_J=[]; LS_A=[]; wLS=[];

% -----------------------------------------------------------------
% Communication optim_fgh with optim_g and optim_H. optim_c w optim_dc
global NLP_xg NLP_g NLP_xdc NLP_dc NLP_xH NLP_H

glbSave(d).NLP_xg      = NLP_xg;
glbSave(d).NLP_g       = NLP_g;
glbSave(d).NLP_xdc     = NLP_xdc;
glbSave(d).NLP_dc      = NLP_dc;
glbSave(d).NLP_xH      = NLP_xH;
glbSave(d).NLP_H       = NLP_H;

% -----------------------------------------------------------------
% Communication optim_fgh with optim_g and optim_H. optim_c w optim_dc
global NLP_d2L NLP_xd2L NLP_lamd2L
glbSave(d).NLP_d2L     = NLP_d2L;
glbSave(d).NLP_xd2L    = NLP_xd2L;
glbSave(d).NLP_lamd2L  = NLP_lamd2L;

% -----------------------------------------------------------------
% Communication qp_f/g 
global QP_x QP_Fx

glbSave(d).QP_x        = QP_x;
glbSave(d).QP_Fx       = QP_Fx;

% -----------------------------------------------------------------
global SEP_z SEP_Jz
glbSave(d).SEP_z        = SEP_z;
glbSave(d).SEP_Jz       = SEP_Jz;

SEP_z=[]; SEP_Jz=[];

% -----------------------------------------------------------------
global JTol gTol 
glbSave(d).JTol  = JTol;
glbSave(d).gTol  = gTol;

% -----------------------------------------------------------------
global n_f n_g n_H n_c n_dc n_d2c n_r n_J n_d2r

glbSave(d).n_f       = n_f;
glbSave(d).n_g       = n_g;
glbSave(d).n_H       = n_H;
glbSave(d).n_c       = n_c;
glbSave(d).n_dc      = n_dc;
glbSave(d).n_d2c     = n_d2c;
glbSave(d).n_r       = n_r;
glbSave(d).n_J       = n_J;
glbSave(d).n_d2r     = n_d2r;

% Init counters to 0
n_f=0; n_g=0; n_H=0; n_c=0; n_dc=0; n_d2c=0; n_r=0; n_J=0; n_d2r=0;

% -----------------------------------------------------------------
% Used by pbuild to build search directions and line search steps from x points
global BUILDP p_dx alphaV X_min X_max F_X X_OLD X_NEW pLen

glbSave(d).BUILDP     = BUILDP;
glbSave(d).p_dx       = p_dx;
glbSave(d).alphaV     = alphaV;
glbSave(d).X_min      = X_min;
glbSave(d).X_max      = X_max;
glbSave(d).F_X        = F_X;
glbSave(d).X_OLD      = X_OLD;
glbSave(d).X_NEW      = X_NEW;
glbSave(d).pLen       = pLen;

BUILDP=[]; p_dx=[]; alphaV=[]; X_min=[]; X_max=[]; X_OLD=[]; X_NEW=[]; pLen=[];

% -----------------------------------------------------------------
% optType defines type of optimization problem
% probType defines problem type, same as optType, or a simpler problem type
% solvType defines solver type

global optType probType

glbSave(d).optType   = optType;
glbSave(d).probType  = probType;

probType=[]; 

% MODIFICATION LOG
%
% 981108   hkh   Written
% 981123   mbk   New globals FD_hRelJ and FD_hRelH
% 990823   hkh   Add QP_x and QP_Fx and optim globals
% 000911   hkh   Another solution for step sizes in derivative estimation
% 000925   hkh   Add optType
% 000116   hkh   Add JTol and gTol used for numerical differences
% 010324   hkh   Add US_A, US_B for user communication
% 020409   hkh   Add NARG for saving of result of xnargin
% 031006   ango  Add FDVar for numerical differences
% 040410   hkh   Remove FDVar, now Prob.FDVar
% 040410   hkh   Remove FD_hRel FD_hRelJ FD_hRelH TIME0 TIME1
% 020418   hkh   Add NARGO for saving of result of nargout
% 050223   frhe  Add d2L globals
% 061212   med   ADMAT removed