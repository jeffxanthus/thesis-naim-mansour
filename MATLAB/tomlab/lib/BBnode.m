% -------------------------------------------------------------------------
function [Node,L] = BBnode(NodeSel,L,pred,Depth,Icomp,fIPMin)
% -------------------------------------------------------------------------
% -------------------------------------------------
% BBnode implements node selection strategies in Branch and Bound algorithms
%
% INPUT:
%   NodeSel   Node selection method in branch and bound
%           = 0 Depth First. Priority on  nodes with more integer components.
%           = 1 Breadth First. Priority on  nodes with more integer components.
%           = 2 Depth First. When integer solution found, use NodeSel = 1 (default)
%           = 3 Pure LIFO (Last in, first out) Depth First
%           = 4 Pure FIFO (First in, first out) Breadth First
%           = 5 Pure LIFO Depth First. When integer solution found, use NodeSel 4
%
%   L         Node list
%   pred      Preceding node in the tree for each node
%   Depth     Depth in the tree
%   Icomp     Number of integer components in each computed node
%   fIPMin    Current best integer feasible solution found
%
% OUTPUT:
%   Node      Selected node in L
%   L         New node list L, without Node

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Oct 15, 2009.   Last modified Oct 16, 2009.

% Problem selection and relaxation
if NodeSel == 0 | (NodeSel == 2 & isinf(fIPMin)) | (NodeSel == 6 & isinf(fIPMin))
   % Find nodes with max depth, deepest down in the tree
   % Select the node with most integer components, first one if ties
   Lev    = Depth(L);
   LevMax = max(Lev);
   ix     = find(Lev==LevMax);
   if length(ix) > 1
      IC  = Icomp(pred(L(ix)));
      [ICmax, ICidx] = max(IC);
      i   = L(ix(ICidx));
   else
      i   = L(ix);
   end

   % fprintf('i %d OLD LIFO i %d\n',i,L(length(L)));
   L = setdiff(L,i);
elseif NodeSel == 1 | NodeSel == 2
   % Find nodes with min depth, highest up in the tree
   % Select the node with most integer components, first one if ties
   Lev    = Depth(L);
   LevMin = min(Lev);
   ix     = find(Lev==LevMin);
   if length(ix) > 1
      IC = Icomp(pred(L(ix)));
      [ICmax, ICidx] = max(IC);
      i  = L(ix(ICidx));
   else
      i  = L(ix);
   end

   % fprintf('i %d OLD FIFO i %d\n',i,L(1));
   L = setdiff(L,i);
elseif NodeSel==3 | (NodeSel==5 & isinf(fIPMin))
   % Old LIFO depth search, last in, first out
   i = L(length(L));
   L = L(1:length(L)-1);
elseif NodeSel==4 | NodeSel==5
   % Old FIFO breadth search, first in, first out
   i = L(1);
   L = L(2:length(L));
elseif NodeSel==6 
   % HKH - Tried new strategy
   IC                 = Icomp(pred(L));
   [ICmax, ICidx]     = max(IC);
   if ICmax == nI-1
      i               = L(ICidx);
   else
      ix = find(IC==ICmax);
      if length(ix) == 1
         i            = L(ICidx);
      else
         [Fmin, fidx] = min(f_min(L));
         i            = L(fidx);
      end
   end
end
Node = i;

% MODIFICATION LOG
%
% 091015  hkh  Written
