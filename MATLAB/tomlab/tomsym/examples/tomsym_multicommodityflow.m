% Minimum Cost Multi-Commodity Flow Problem (tomSym example) 
%
% The multi-commodity flow problem is a network flow problem with multiple
% commodities (or goods) flowing through the network.
%
% In this variation of the problem, the cost of flow through the network is
% minimized, while meeting a constraint on the demand at the output node.
%
% Reference: http://en.wikipedia.org/wiki/Multi-commodity_flow_problem

% Network description:
% Each row of the "edges" matrix corresponds to one edge in the network.
% The first row is the input node and the last row is the output node.
% (In a real problem, there could be tens of thousands of edges.)
edges = [...
    1 2 10 1 1 1 1.7
    1 3  3 0 1 0 0.7
    1 3  3 0 0 1 2.0
    1 4  5 1 1 1 9.9
    2 3  8 1 1 0 1.2
    3 4 10 1 1 1 1.0
    ];

from = edges(:,1);  % Start node for each edge.
to = edges(:,2);  % End node for each edge.
hasflow = logical(edges(:,4:6)); % The commodity is allowed to flow 
                                 % on the egge (1) or not (0).
ub = edges(:,3); % Upper bound on total flow through the edge.
cst = edges(:,7); % Cost per unit for flow through the edge.

demand = [3 4 5];

Nc = size(hasflow,2); % Number of commodities
Ne = size(hasflow,1); % Number of edges
Nn = max(to(:));  % Number of nodes

% Set up a symbolic decsion variable for each commodity on each edge where
% it is allowed to flow.
flowv = tom('flowv',sum(hasflow(:)),1,'integer');

% Make a sparse symbolic matrix, filling in with zeros for
% edges/commodities that are disallowed.
[i,j] = find(hasflow);
flow = sparse(i,j,flowv,Ne,Nc);

% Non-negativity constraint
c1 = ( flowv >= 0 );

% Capacity constraints
flowsum = sum(flow,2); % Sum of all commodities through each edge.
c2 = ( flowsum <= ub );

% Flow conservation constraints 
Ai = sparse(from,(1:Ne)',ones(Ne,1),Nn,Ne);
Ao = sparse(to,(1:Ne)',ones(Ne,1),Nn,Ne);
A = Ai-Ao;
A = A(2:end-1,:); % First and last nodes are sink and source, 
                  % so they have no conservation constraint.
c3 = cell(1,Nc);
for i=1:Nc
    c3{i} = ( A*flow(:,i) == 0 );
end

% Demand constraints
c4 = cell(1,Nc);
for i=1:Nc
    c4{i} = ( Ao(end,:)*flow(:,i) == demand(i) );
end

% Objective function
total_cost = sum(cst.*flowsum);

% Starting guess
guess = ( flowv == 0 );

% Solve the problem
options = struct;
options.name = 'multi-commodity mincost';
options.solver = 'cplex'; % gurobi also works

solution = ezsolve(total_cost, {c1,c2,c3,c4}, guess, options);

% Substitute the solution into our flow variable, to get a numeric result.
disp('Flows per edge and commodity:');
disp(full(subs(flow,solution)));
