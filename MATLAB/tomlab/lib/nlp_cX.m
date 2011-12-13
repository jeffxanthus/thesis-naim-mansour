% nlp_cX.m
%
% function [c, Pdim]=nlp_cX(x, Prob, varargin)
%
% Calls the interface routine nlp_cF.m
%
% nlp_cX also explicitely add linear constraints
%
% In Pdim the different constraint dimensions are stored
%
% (1) Total number of constraints (not bounds) (=mcon)
% (2) Number of linear and nonlinear equalities (=me)
% (3) Number of linear equalities (=meq)
% (4) Number of linear inequalities (=mleq)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Aug 25, 1998. Last modified Oct 19, 2000.

function [c, Pdim]=nlp_cX(x, Prob, varargin)
[c, Pdim]= nlp_cF( x, Prob, varargin{:});