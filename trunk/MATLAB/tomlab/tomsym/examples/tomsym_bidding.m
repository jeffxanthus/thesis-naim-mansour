%% Bid Evaluation
% TomSym implementation of GAMS Example (BID,SEQ=19)
%
% A company obtains a number of bids from vendors for the supply
% of a specified number of units of an item. Most of the submitted
% bids have prices that depend on the volume of business.
%
% Bracken, J, and McCormick, G P, Chapter 3. In Selected Applications of
% Nonlinear Programming. John Wiley and Sons, New York, 1968, pp. 28-36.
%
% v: vendors (a, b, c, d, e)
%
% s: segments (1-5)
%
% vs(v,s): vendor bit possibilities
%
% cl: column labels (setup, price, q-min, q-max)

% scalar req requirements
req = 239600.48;

% Matrix bid(v,s,cl), with bid data
bid = zeros(5,5,4);
bid(1,1,:) = [3855.84     61.150     0        33000];
bid(2,1,:) = [125804.84   68.099     22000    70000];
bid(2,2,:) = [0           66.049     70000    100000];
bid(2,3,:) = [0           64.099     100000   150000];
bid(2,4,:) = [0           62.119     150000   160000];
bid(3,1,:) = [13456.00    62.190     0        165600];
bid(4,1,:) = [6583.98     72.488     0        12000];
bid(5,1,:) = [0           70.150     0        42000];
bid(5,2,:) = [0           68.150     42000    77000];

% Get minimum domains and ripple total cost up the segments
vs = bid(:,:,4);

% bid(v,s+1,'setup') = bid(v,s,'setup') +
% bid(v,s,'q-max')*(bid(v,s,'price')-bid(v,s+1,'price')))
for i=1:5
    for j=2:5
        bid(i,j,1) = bid(i,j-1,1)+bid(i,j-1,4)*(bid(i,j-1,2)-bid(i,j,2));
    end
end

% Variables

% pl(v,s): Purchase level
toms 5x5 pl

% plb(v,s): Purchase decision (binary)
toms 5x5 integer plb

% Bounds on binary (integer variables)
cbnd = {0 <= plb <= 1};

% Index for variables that are known to be zero
idx = vs==0;

% Set bounds to force zero
cbnd = {cbnd; pl(vs==0) == 0; plb(vs==0) == 0};

% Demand constraint
eq1 = {req == sum(sum(pl))};

% Cost definition
bid1 = bid(:,:,1);
bid1(idx) = 0;
bid2 = bid(:,:,2);
bid2(idx) = 0;
cost = sum(sum(bid2.*pl + bid1.*plb));

% Min purchase (pl(vs) >= bid(vs,'q-min')*plb(vs))
bid3 = bid(:,:,3);
bid3(idx) = 0;
eq2 = {pl >= bid3.*plb};

% Max purchase (pl(vs) <= bid(vs,'q-max')*plb(vs))
bid4 = bid(:,:,4);
bid4(idx) = 0;
eq3 = {pl <= bid4.*plb};

% At most one deal
eq4 = {sum(plb,2) <= 1};

% Constraints
cons = {cbnd,eq1,eq2,eq3,eq4};

solution = ezsolve(cost,cons);