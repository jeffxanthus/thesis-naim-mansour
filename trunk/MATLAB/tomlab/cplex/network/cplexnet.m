% TOMLAB /CPLEX Network Solver
%
% -----------------------------------------------------------------------
%
%  cplexnet.m solves network planning problems:
%
%  Minimize (or maximize) SUM (obj.*x) for all arcs.
%
%  Subject to:
%  For all nodes in the network:
%  SUM(x) for all arcs in tail - SUM(x) for all arcs in heads = Supply.
%  lb <= x <= ub for all arcss
% 
%  That is, for each node, the net flow entering and leaving the node must equal 
%  its supply value, and all flow values must be within their bounds. The solution 
%  of a network-flow problem is an assignment of flow values to arcs (that is, 
%  the modeling variables) to satisfy the problem formulation. A flow that satisfies 
%  the constraints and bounds is feasible. 
%  
%  Finds minimal-cost flow through a network. The network has N nodes and a set of arcs
%  connecting the nodes. An arc a in the set A is an ordered pair (i, j) where i and j 
%  are nodes in the set N; node i is called the tail or the from-node and node j is called 
%  the head or the to-node of the arc a. Not all the pairs of nodes in a set N are 
%  necessarily connected by arcs in the set A. More than one arc may connect a pair of 
%  nodes; in other words, a1 = (i, j) and a2 = (i, j) may be two different arcs in A, 
%  both connecting the nodes i and j in N. 
%
%  Each arc may be associated with four values:
%
%  x   - The flow passing the arc form its tail to its head. Neg. values
%        allowed.
%  lb  - Minimum flow through the arc.
%  ub  - Maximum flow thorugh the arc.
%  obj - The objective, contribution to objective function for one unit
%        flow.
%
%  Each node is associated with one value:
%
%  supply - The supply value at the node. A node with a positive value is
%  called a supply node (source), if negative a demand node or sink. If the
%  value is 0, a transshipment node. The sum of all supplies must match the
%  sum of all demands, or the problem is infeasible.
% 
%  tail is the set of arcs whose tails are node n. head is the set of arcs
%  whose heads are node n. 
%
% ------------------------------------------------------------------------
%
% function  [x, slack, v, rc, f_k, Inform, Iter] = cplexnet(obj, ub, lb, tail, head,
%           supply, callback, PriLev, BIG, cpxControl, logfile, savefile, savemode, netfile);
%
% INPUT:  
% obj       Objective function cost coeffs. One for each arc.
% ub        Upper bounds for each arc.
% lb        Lower bounds for each arc.
% 
% tail      Tail of the arc (start node)
% head      Head of the arc (end node)
% supply    Supply and demand nodes, + supply, - demand.
%
%    The following parameters are optional:
%
% callback  Logical scalar defining if callback is used in CPLEX
%  callback = 1 activates the callback. See TOMLAB /CPLEX User's Guide
%  The callback calls the m-file specified below. The user may edit this file,
%  or make a new copy, which is put before in the Matlab path.
%
% callback(1)  cpxcb_NET.m        Network simplex callback
%
% PriLev    Printing level in the CPLEX m-file and CPLEX C-interface.
%           = 0    Silent
%           = 1    Warnings and Errors
%           = 2    Summary information
%           = 3    More detailed information
%
%           > 10   Pause statements, and maximal printing (debug mode)
%
% BIG       Defines default lower and upper bounds, default 1E20.
%
% cpxControl Structure, where the fields are set to the CPLEX 
%            parameters that the user wants to specify values for.
%            The following parameters are the only ones of general
%            interest. Default values are recommended:
%
%  cpxControl.NETITLIM: Limits the number of iterations that the network 
%             optimizer performs. Default BIGINT.
% 
%  cpxControl.NETEPOPT: Optimality tolerance for the network optimizer. 
%             The optimality tolerance specifies the amount a reduced cost 
%             may violate the criterion for an optimal solution.
%             Default 1e-6. Valid values from 1e-11 to 1e-1.             
%
%  cpxControl.NETEPRHS: Feasibility tolerance for the network optimizer. 
%             The feasibility tolerance specifies the degree to which a problem's
%             flow value may violate its bounds. This tolerance influences the
%             selection of an optimal basis and can be reset to a higher value 
%             when a problem is having difficulty maintaining feasibility
%             during optimization. You may also wish to lower this tolerance 
%             after finding an optimal solution if there is any doubt that 
%             the solution is truly optimal. If the feasibility tolerance is 
%             set too low, CPLEX may falsely conclude that a problem is infeasible.
%             If you encounter reports of infeasibility in the optimization, 
%             a small adjustment in the feasibility tolerance may improve performance.
%             Default 1e-6. Valid values from 1e-11 to 1e-1.
% 
%  cpxControl.NETPPRIIND: Pricing algorithm for the network optimizer. 
%             On the rare occasions when the network optimizer seems to 
%             take too long to find a solution, you may want to change 
%             the pricing algorithm to try to speed up computation. All 
%             the choices use variations of partial reduced-cost pricing.
%
%             NETPPRIIND = 0: automatic, default (same as 3)
%             NETPPRIIND = 1: Partial pricing.
%             NETPPRIIND = 2: Multiple partial pricing.
%             NETPPRIIND = 3: Multiple partial pricing with sorting.
% 
%  cpxControl.NETFIND: The CPLEX network extractor searches an LP constraint 
%             matrix for a submatrix with the following characteristics: 
%             - the coefficients of the submatrix are all 0 (zero), 1 (one), 
%               or -1 (minus one); 
%             - each variable appears in at most two rows with at most one 
%               coefficient of +1 and at most one coefficient of -1. 
%             CPLEX can perform different levels of extraction. The level 
%             it performs depends on the NETFIND parameter. 
%
%             NETFIND = 1: CPLEX extracts only the obvious network; it 
%             uses no scaling; it scans rows in their natural order; it 
%             stops extraction as soon as no more rows can be added to the 
%             network found so far. 
%             NETFIND = 2: Default. CPLEX also uses reflection scaling (that is, it 
%             multiplies rows by -1) in an attempt to extract a larger
%             network.
%             NETFIND = 3: CPLEX uses general scaling, rescaling both rows and 
%             columns, in an attempt to extract a larger network. 
%
%             In terms of total solution time expended, it may or may not be 
%             advantageous to extract the largest possible network. Characteristics 
%             of your problem will determine the tradeoff between network size and 
%             the number of simplex iterations required to finish solving the model 
%             after solving the embedded network. 
% 
%             Even if your problem does not conform precisely to network conventions, 
%             the network optimizer may still be advantageous to use. When it is possible 
%             to transform the original statement of a linear program into network 
%             conventions by these algebraic operations: 
%
%             - changing the signs of coefficients.
%             - multiplying constraints by constants.
%             - rescaling columns.
%             - adding or eliminating redundant relations.
%
%             then CPLEX will carry out such transformations automatically if 
%             you set the NETFIND parameter appropriately. 
%
%  cpxControl.PREPASS: If your LP problem includes network structures, there is a 
%             possibility that CPLEX preprocessing may eliminate those structures 
%             from your model. For that reason, you should consider turning off 
%             preprocessing before you invoke the network optimizer on a problem. 
%
%             PREPASS = -1: Default. Determined automatically.
%             PREPASS = 0: Do not use Presolve
%
% logfile   Name of file to write CPLEX log to. If empty, no log is
%           written. 
%
% savefile  Name of file to write CPLEX problem just prior to calling the 
%           CPLEX solver. If empty, nothing is written. Also see the savemode 
%           parameter below.
%
% savemode  Integer flag indicating which format to use for the save file.
%           The following values are possible:
%
%        1: SAV  Binary SAV file
%        2: MPS  MPS format, original format
%        3: LP   CPLEX LP format, original format
%        4: RMP  MPS file with generic names
%        5: REW  MPS file with generic names
%        6: RLP  LP  file with generic names
%
%            The SAV format is a binary format suitable for submission to
%            ILOG help desk.
%
% netfile   MATLAB mat file defining the problems first 6 inputs.
%
%           For example netfile = 'nettest1.mat'.
%
% ------------------------------------------------------------------------------
%
% OUTPUT: 
%
% x         Solution vector with decision variable values (n x 1 vector)
% slack     Slack variables (m x 1 vector)
% v         Lagrangian multipliers (dual solution vector) (m x 1 vector)
% rc        Reduced costs. Lagrangian multipliers for simple bounds on x.
% f_k       Objective function value at optimum
%
% Inform    Result of CPLEX run: 
%   
%       1 Optimal solution is available
%       2 Model has an unbounded ray
%       3 Model has been proved infeasible
%       4 Model has been proved either infeasible or unbounded
%       5 Optimal solution is available, but with infeasibilities after unscaling
%       6 Solution is available, but not proved optimal, due to numeric difficulties 
%      10 Stopped due to limit on number of iterations
%      11 Stopped due to a time limit
%      12 Stopped due to an objective limit
%      13 Stopped due to a request from the user
%      
%      32 Converged, dual feasible, primal infeasible
%      33 Converged, primal feasible, dual infeasible
%      34 Converged, primal and dual infeasible
%      35 Primal objective limit reached
%      36 Dual objective limit reached
%      37 Primal has unbounded optimal face
%      38 Non-optimal solution found, primal-dual feasible
%      39 Non-optimal solution found, primal infeasible
%      40 Non-optimal solution found, dual infeasible
%      41 Non-optimal solution found, primal-dual infeasible
%      42 Non-optimal solution found, numerical difficulties
%      43 Barrier found inconsistent constraints
%
%     101 Optimal integer solution found
%     102 Optimal sol. within epgap or epagap tolerance found
%     103 Solution is integer infeasible
%     104 The limit on mixed integer solutions has been reached 
%     105 Node limit exceeded, integer solution exists
%     106 Node limit exceeded, no integer solution
%     107 Time limit exceeded, integer solution exists
%     108 Time limit exceeded, no integer solution
%     109 Terminated because of an error, but integer solution exists.
%     110 Terminated because of an error, no integer solution 
%     111 Limit on tree memory has been reached, but an integer solution exists 
%     112 Limit on tree memory has been reached; no integer solution 
%     113 Stopped, but an integer solution exists.
%     114 Stopped; no integer solution.
%     115 Problem is optimal with unscaled infeasibilities 
%     116 Out of memory, no tree available, integer solution exists 
%     117 Out of memory, no tree available, no integer solution 
%     118 Model has an unbounded ray 
%     119 Model has been proved either infeasible or unbounded 
%
% Iter    Number of iterations

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2004 by Tomlab Optimization Inc., $Release: 9.0.2$
% Written Apr 16, 2004.    Last modified Mar 22, 2005.

function [x, slack, v, rc, f_k, Inform, Iter] = cplexnet(obj, ub, lb, tail, head, ...
          supply, callback, PriLev, BIG, cpxControl, logfile, savefile, savemode, netfile)
      
if nargin<14
   netfile = [];
     if nargin < 13
        savemode = [];
           if nargin < 12
              savefile = '';
                 if nargin < 11
                    logfile = '';
                      if nargin < 10
                         cpxControl = [];
                           if nargin < 9
                              BIG = [];
                                if nargin < 8
                                   PriLev = [];
                                     if nargin < 7
                                        callback = [];
                                          if nargin < 6
                                             error('cplexnet needs at least 6 arguments');
end, end, end, end, end, end, end, end, end

if isempty(netfile)
    n = max( [length(ub),length(lb),length(obj)] );
    if n==0, error('cplexnet: cannot determine problem dimension'); end
    if sum(supply)~=0, error('cplexnet: supply vector must sum to zero'); end
    if min(tail) < 1, error('cplexnet: tail has node indices below 1'); end
    if min(head) < 1, error('cplexnet: head has node indices below 1'); end
    if max(tail) > length(supply), error('cplexnet: tail has indices out of range'); end
    if max(head) > length(supply), error('cplexnet: head has indices out of range'); end
end

% Check on settings:
if ~isempty(cpxControl)
    if isfield(cpxControl,'NETITLIM')
        if cpxControl.NETITLIM < 0, error('cplexnet: NETITLIM less than zero'); end
    end
    if isfield(cpxControl,'NETEPOPT')
        if cpxControl.NETEPOPT < 1e-11, error('cplexnet: NETEPOPT too low'); end
        if cpxControl.NETEPOPT > 1e-1, error('cplexnet: NETEPOPT too high'); end
    end
    if isfield(cpxControl,'NETEPRHS')
        if cpxControl.NETEPRHS < 1e-11, error('cplexnet: NETEPRHS too low'); end
        if cpxControl.NETEPRHS > 1e-1, error('cplexnet: NETEPRHS too high'); end
    end
    if isfield(cpxControl,'NETPPRIIND')
        if ~any(cpxControl.NETPPRIIND == [0 1 2 3]), error('cplexnet: NETPPRIIND invalid'); end
    end
    if isfield(cpxControl,'NETFIND')
        if ~any(cpxControl.NETFIND == [0 1 2 3]), error('cplexnet: NETFIND invalid'); end
    end
end

obj = obj(:);
ub  = ub(:);
lb  = lb(:);
tail = tail(:); % Indices moved
head = head(:); % Indices moved
supply = supply(:);

if isempty(PriLev), PriLev = 0; end
if isempty(callback), callback=0; end
if isempty(BIG), BIG=1e20; end

% All other inputs checked by nargin.

try
   [x, slack, v, rc, f_k, Inform, Iter] = cplexnetmex(obj, ub, lb, tail, head,...
     supply, callback, PriLev, BIG, cpxControl, logfile, savefile, savemode, netfile);
catch
   V = version;
   if(strcmpi(V(1:3),'6.0') | strcmpi(V(1:3),'6.1'))
      l = lasterr;
      if(strcmp(l,'Invalid MEX-file'))
         tomlabsharederror;
      else
         error(l);
      end
   else
      l=lasterror;
      if(strcmp(l.identifier,'MATLAB:invalidMEXFile'))
         tomlabsharederror;
      else
         rethrow(l);
      end
   end
end

% MODIFICATION LOG:
%
% 040416 med  Created
% 040502 med  Netfile check added for errors
% 050117 med  Revised with mlint
