% probInit.m :
%
% General initialization routine for TOMLAB optimization problems,
% when a file in the TOMLAB Init File format is used
%
% function [Prob] = probInit (probFile, P, ask, Prob);
%
% INPUT:
% probFile Name of problem definition Init File, as a string, e.g. 'con_prob'.
% P        Problem number
%          If P=0, call Init File to get the defined problems and let the user
%          choose a problem by displaying a menu
%          if isempty(P), just return the number of problems
% ask      1:  ask questions; 0: use defaults; -1: use values in Prob.uP
%          11: ask questions in the GUI window
%
%          If ask >= 2, ask = dimension of problem, for init files with 
%          variable problem dimension, and n as a 3rd input parameter
%          example: Prob = probInit('glbv_prob',1,50);
%          Problem 1 in glbb_prob is defined with dimension 50
%
%          If isempty(ask) then [if isempty(Prob.uP),ask=1, else ask=-1];
%
% Prob     Either problem structure with all parameters defining the problem
%          Or the numerical vector uP, to be set as Prob.uP and used in the
%          definition of the problem in the Init File
%
% OUTPUT:
% Prob     Problem structure or
%          Number of problems, if isempty(P)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.2.0$
% Written June 16, 1999. Last modified June 24, 2008.

function Prob = probInit (probFile, P, ask, Prob)

if nargin < 4
   Prob=[];
   if nargin < 3
      ask=[];
      if nargin < 2
         P=[];
         if nargin < 1
            probFile='con_prob';
         end
      end
   end
end

if isempty(P)   % Only return the number of predefined problems
   probList=feval(probFile);
   nProbs=size(probList,1);
   Prob=nProbs;
   return
end

if isempty(ask) 
   if isstruct(Prob)
      if isempty(Prob.uP) 
         ask=1; 
      else
         ask=-1; 
      end
   else
      ask=-1; 
   end
end

if P(1)<=0  % Give menu to choose problem from
   P=strmenu('Choice of test function ',feval(probFile,[]));
elseif ask(1) > 1
   % Do nothing
elseif ~isstruct(Prob) & ~isempty(Prob)
   % Assume uP is sent instead of Prob
   uP            = Prob;
   probList      = feval(probFile);
   Prob          = ProbDef;
   Prob.P        = P; 
   Prob.probFile = []; 
   Prob.Name     = probList(P,:); 
   Prob.uP       = uP;
end


[probList, Prob] = feval(probFile, P, ask, Prob);

global probType

probType=Prob.probType;

Prob.probFile=probFile;

% MODIFICATION LOG:
%
% 980825  hkh  Changed error in comments. Changed call to usr_prob-routine.
% 980910  mbk  Do not set Prob to [] on line 72.
% 981026  hkh  Changed ProbFile to probFile. Set field Prob.probFile to probFile
% 990617  hkh  Change comments
% 990622  hkh  Add first argument (empty) to initial call to probFile,usr_prob
% 001012  hkh  Delete usr_prob possibility, not needed any longer
% 050228  hkh  Change prob.uP to Prob.uP
% 050302  hkh  Prob or uP as 4th input, to enable Prob.uP to be set
% 080624  hkh  Ask optionally the problem dimension, used e.g. by glbv_prob
