% function Prob = refinery(demand, supply, capacity, costs, supplycomp, numvbls, octane,...
%   vappres, volatility, sulfur, reformer, cracker, petrolspec_L, petrolspec_U, dieselspec_L,...
%   dieselspec_U);
%
% Creates a TOMLAB LP problem for refinery production
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob         A TOMLAB problem defining the problem (type LP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.1.0$
% Written Oct 7, 2005.   Last modified Jan 4, 2006.

function Prob = refinery(demand, supply, capacity, costs, supplycomp, numvbls, octane,...
   vappres, volatility, sulfur, reformer, cracker, petrolspec_L, petrolspec_U, dieselspec_L,...
   dieselspec_U)

if nargin < 15
   error('The function requires 13 inputs');
end

demand       = demand(:);
supply       = supply(:);
capacity     = capacity(:);
costs        = costs(:);
octane       = octane(:);
vappres      = vappres(:);
volatility   = volatility(:);
sulfur       = sulfur(:);
reformer     = reformer(:);
cracker      = cracker(:);
petrolspec_L = petrolspec_L(:);
petrolspec_U = petrolspec_U(:);
dieselspec_L = dieselspec_L(:);
dieselspec_U = dieselspec_U(:);

% PRODUCTS
% butane, petrol, diesel, heating        (FINAL PRODUCTS, 4) 4
% petbutane, reformerate, petcrknaphtha  (IPETROL, 3) 7
% dslgasoil, dslcrknaphtha, dslcrkgasoil (IDIESEL, 3) 10
% hogasoil, hocrknaphtha, hocrkgasoil    (IHO, 3) 13
% distbutane, naphtha, residue, gasoil   (IDIST, 4) 17
% refbutane, reformate                   (IREF, 2) 19
% crknaphtha, crkgasoil                  (ICRACK, 2) 21
% crude1, crude2                         (CRUDES, 2) 23

n1 = length(demand);
n3 = length(supply);
n2 = numvbls-n1-n3;
n  = numvbls;

% FORMULATE PROBLEM

% Objective function
c  =  [zeros(n1+n2,1);costs(1:n3,1)]; % Cost of crudes
c(15:17,1) = costs(n3+1:n3+3);        % the 3 final ones need further processing

% Production constraint
idist = 4;
A1 = zeros(idist,n); % 4 IDIST, starts at 14, crudes are 22 and 23
for i=1:idist
   A1(i,[13+i, 22, 23]) = [1,-supplycomp(:,i)'];
end
b_L1 = -inf*ones(idist,1);
b_U1 = zeros(idist,1);

% Reformer constraint
iref = 2;
A2   = zeros(iref,n);
for i=1:iref
   A2(i,[17+i, 15]) = [1, -reformer(i)]; %Naphtha
end
b_L2 = -inf*ones(iref,1);
b_U2 = zeros(iref,1);

% Cracker constraint
icrack = 2;
A3   = zeros(icrack,n);
for i=1:icrack
   A3(i,[19+i, 16]) = [1, -cracker(i)]; %Residue
end
b_L3 = -inf*ones(icrack,1);
b_U3 = zeros(icrack,1);

% Product equality constraints
A4 = zeros(3,n);
A4(1, [20, 7, 12, 9]) = [1 -1 -1 -1];
A4(2, [21, 13, 10])   = [1 -1 -1];
A4(3, [17, 11, 8])    = [1 -1 -1];
b_L4 = zeros(3,1);
b_U4 = inf*ones(3,1);

% Final produts equalities
A5 = zeros(4,n);
A5(1,[1, 14, 18, 5]) = [1 -1 -1 1];
A5(2,[2, 5:7])       = [1 -1 -1 -1];
A5(3,[3, 8:10])      = [1 -1 -1 -1];
A5(4,[4, 11:13])     = [1 -1 -1 -1];
b_L5 = zeros(4,1);
b_U5 = b_L5;

% Octane constraint
A6 = zeros(4,n);
A6(1,[5:7, 2]) = [octane', -petrolspec_L(1)];
A6(2,[5:7, 2]) = [vappres', -petrolspec_U(1)];
A6(3,[5:7, 2]) = [volatility', -petrolspec_L(2)];
A6(4,[8:10, 3]) = [sulfur', -dieselspec_U(1)];
b_L6 = [0;-inf;0;-inf];
b_U6 = [inf;0;inf;0];

% Reformerate are the same
A7 = zeros(1,n);
A7(1, [6, 19]) = [1 -1];
b_L7 = 0;
b_U7 = 0;

% Merge constraints
A = [A1;A2;A3;A4;A5;A6;A7];
b_L = [b_L1;b_L2;b_L3;b_L4;b_L5;b_L6;b_L7];
b_U = [b_U1;b_U2;b_U3;b_U4;b_U5;b_U6;b_U7];

x_L = zeros(n,1);
x_U = inf*ones(n,1);
x_U([15, 16, 17],1) = capacity; %Not same order as orig model
x_U([22, 23],1)     = supply;   %Available crude
x_L([1:4], 1)       = demand;   %Demands

Prob = lpAssign(c, A, b_L, b_U, x_L, x_U, [], 'Refinery');

% MODIFICATION LOG
%
% 051007 med   Created
% 060104 med   Updated and corrected