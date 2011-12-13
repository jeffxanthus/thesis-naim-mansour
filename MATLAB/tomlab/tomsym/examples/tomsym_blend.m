%% Blending Problem I
% TomSym implementation of GAMS Example (BLEND,SEQ=2)
%
% A company wishes to produce a lead-zinc-tin alloy at minimal cost.
% The problem is to blend a new alloy from other purchased alloys.
%
% Dantzig, G B, Chapter 3.4. In Linear Programming and Extensions.
% Princeton University Press, Princeton, New Jersey, 1963.
%
% alloy: products on the market (a, b, c, d, e, f, g, h, i)
%
% elem: required elements (lead, zinc, tin)

% Composition data (percent)
compdat = [10   10   40   60   30   30   30   50   20;
    10   30   50   30   30   40   20   40   30;
    80   60   10   10   40   30   50   10   50];

% Alloy price
price   = [4.1  4.3  5.8  6.0  7.6  7.5  7.3  6.9  7.3];

% Parameters rb(elem)  required blend / lead = 30, zinc = 30, tin = 40 /
%            ce(alloy) composition error (pct-100);

% Required blend (lead = 30, zinc = 30, tin = 40)
rb = [30;30;40];

% Variables

% Purchase of alloy (pounds)
toms 9x1 alloy

% Variable has to be positive
bnd = {alloy >= 0};

% Purchase constraint
pc = {rb == compdat*alloy};

% Material balance
mb = {1 == sum(alloy)};

% Accounting: total cost
phi = price*alloy;

options = struct;
options.name = 'Blend';
solution = ezsolve(phi,{pc, mb, bnd},{},options);