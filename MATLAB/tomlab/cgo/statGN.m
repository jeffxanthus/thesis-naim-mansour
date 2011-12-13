% ====================================================================
% function [onB, doX, doM] = statGN(x_L,x_U,xLoc,xMin,X,epsX,...
%    IntVars, Reals)
% ====================================================================

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2008 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written Nov 24, 2006. Last modified Feb 24, 2008.

function [onB, doX, doM] = statGN(x_L,x_U,xLoc,xMin,X,epsX,...
   IntVars, Reals)

if nargin < 8
   if nargin < 7
      IntVars = [];
   end
   if ~isempty(IntVars)
      IV          = ones(length(x_L),1);
      IV(IntVars) = 0;
      Reals       = find(IV);
   end
end

if isempty(IntVars)
   onB = nOnBound(xLoc,x_L,x_U,epsX);
   doX = min(tomsol(30,xLoc,X));
   doM = min(tomsol(30,xLoc,xMin));
else
   onB = nOnBound(xLoc(Reals),x_L(Reals),x_U(Reals),epsX);
   iv = double(xLoc(IntVars(1)) ~= X(IntVars(1),:));
   for i = 2:length(IntVars)
       iv = iv + double(xLoc(IntVars(i)) ~= X(IntVars(i),:));
   end
   doX(2) = min(iv);
   doM(2) = sum(xLoc(IntVars)~=xMin(IntVars));
   
   if(isempty(Reals))
       doM(1) = 0;
       doX(1) = 0;
   else
       doM(1) = min(tomsol(30,xLoc(Reals),xMin(Reals)));
       %ix     = find(iv == doX(2));
       %doX(1) = min(tomsol(30,xLoc(Reals),X(Reals,ix)));
       doX(1) = min(tomsol(30,xLoc(Reals),X(Reals, iv == doX(2))));
   end

end

