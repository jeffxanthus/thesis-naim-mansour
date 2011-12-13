% TOMLAB - MAD Dummy Interface
% Version 7.2 (R7.2.0) 14-Aug-2009
%
% tomlab\mad Dummy routines for MAD TB, used in TOMLAB
%
% Contents     This file
%
% To use automatic differentiation with the MAD TB:
%
% Install MAD TB.
% Follow the instructions in the manual to complete the installation.
%
% Run startupmad before using MAD in TOMLAB.
% Running startupmad will place paths to the MAD routines before
% the path to the dummy routines below.
%
% In Tomlab the flag Prob.ADObj=1 and Prob.ADCons=1 makes Tomlab use MAD.
% MAD supports one level of differentiation. The following settings are
% available:
%
% AD for gradient, nlp_g set Prob.ADObj   =  1
% AD for Hessian,  nlp_H set Prob.ADObj   = -1
%
% AD for Jacobian,       nlp_dc  set Prob.ADCons =  1
% AD for 2nd Cons deriv, nlp_d2c set Prob.ADCons = -1
%
% This flag could be set in the menues, GUI, or directly before by the
% user in the Prob structure before call to a solver or driver routine.
%
% fmad          Creates a MAD object which includes the derivatives