% Check if user parameters uP is properly set for the actual problem
%
% function uP=checkuP(Name, Prob);
%
% Name    Name of test problem to be using the       Starting values
%
% Prob    Structure with problem definition. Used fields:
%         uP     User parameters
%         uPName Problem name set when defining user parameters
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Oct 6, 1999.   Last modified June 22, 1999.

function uP=checkuP(Name, Prob)

% User parameters uP, may be defined globally.

if isfield(Prob,'uP')  
   uP=Prob.uP;
   if isfield(Prob,'uPName')
      if ~isempty(Prob.uPName)
         uPName=deblank(Prob.uPName);
         if ~strcmpi(uPName,deblank(Name)), uP=[]; end
      end
   end
else
   uP=[];
end

% MODIFICATION LOG:
%
% 050302  hkh  Change strcmp to strcmpi, not forcing exakt case match
% 050302  hkh  deblank both uPName and Name, traling blanks not significant