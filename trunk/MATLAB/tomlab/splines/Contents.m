% TOMLAB - Splines Dummy Interface
% Version 7.8 (R7.8.0) 2-Dec-2011 
%
% tomlab\splines  Spline dummy routines for TOMLAB 
%
% If the user has not the SPLINES TB installed, these dummy routines make
% it possible to avoid feval, and compile with MCC.
%
% If the user installs the SPLINES TB, the path to this directory must come
% after in the MATLAB PATH
% 
% In Tomlab these routines are used for three alternative methods to 
% compute numerical differences to approximate first and second derivatives
%
% contents     This file
%
% csapi   Cubic spline interpolant
% csaps   Cubic smoothing spline
% fnder   Returns the first derivative of a spline representation
% fnval   Evaluates a spline function
% spaps   Smoothing spline
