%% Simple Farm Level Model
% TomSym implementation of GAMS Example (DEMO1,SEQ=91)
%
% This version has 7 principal crops and 2 basic inputs, land and labor,
% which are specified on a monthly basis.
%
% Kutcher, G P, Meeraus, A, and O'Mara, G T, Agriculture Sector and Policy
% Models. The World Bank, 1988.
%
% c: crops (wheat, clover, beans, onions, cotton, maize,  tomato)
%
% t: period (jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec)

% landreq(t,c): months of land occupation by crop (hectares)
landreq = [1       1      1      1     0     0     0;
    1      1       1      1      0     0     0;
    1      0.5     1      1      0.5   0     0;
    1      0       1      1      1     0     0;
    1      0       0      0.25   1     0.25  0;
    0      0       0      0      1     1     0;
    0      0       0      0      1     1     0.75;
    0      0       0      0      1     1     1;
    0      0       0      0      1     1     1;
    0      0       0      0      1     0.5   1;
    0.5    0.25    0.25   0.5    0.75  0     0.75;
    1      1       1      1      0     0     0];

% laborreq(t,c) crop labor requirements (man-days per hectare)
laborreq = [1.72   4.5     0.75   5.16    0     0      0;
    0.5    1       0.75   5      0     0      0;
    1      8       0.75   5      5     0      0;
    1      0       16    19.58   5     0      0;
    17.16   0       0    2.42    9    4.3     0;
    2.34   0        0    0       2    5.04    0;
    0      0        0    0       1.5   7.16   17;
    0      0        0    0       2    7.97   15;
    0      0        0    0       1    4.41   12;
    0      0        0    0      26    1.12    7;
    2.43   2.5    7.5   11.16   12    0       6;
    1.35   7.5    0.75   4.68    0     0      0];

% Crop yield (tons per hectare)
yield = [1.5 6.5 1 3 1.5 2 3]';

% Crop prices (dollars per ton)
price = [100 0 200 125 350 70 120]';

% Misc cash costs (dollars per hectare)
misccost = [10 0 5 50 80 5 50]';

% Farm size (hectares)
land   = 4;

% Family labor available (days per month)
famlab = 25;

% Hire-out  wage rate (dollars per day)
owage  = 3;

% Temporary labor wage  (dollars per day)
twage  = 4;

% Number of working days per month
dpm    = 25;

% Cropping activity
toms 7x1 xcrop

% Family labor use, hiring out and temporary labor
toms 12x1 flab fout tlab

% All variable are positive
boxcon = {0 <= xcrop
    0 <= flab
    0 <= fout
    0 <= tlab};

% Land balance
t = size(landreq,1);
eq1 = {};
for i=1:t
    eq1 = {eq1; sum(xcrop.*landreq(i,:)') <= land};
end

% Labor balance
eq2 = {};
for i=1:t
    eq2 = {eq2; sum(xcrop.*laborreq(i,:)') <= flab(i)+tlab(i)};
end

% Family labor balance
eq3 = {};
for i=1:t
    eq3 = {eq3; famlab == flab(i)+fout(i)};
end

% Revenue accounting
revenue = sum(xcrop.*yield.*price);

% Cash cost accounting
mcost = sum(xcrop.*misccost);

% Labor cost accounting
labcost = sum(tlab.*twage);

% Labor income accounting
labearn = sum(fout.*owage);

% Income definition
yfarm = revenue + labearn - labcost - mcost;

solution = ezsolve(-yfarm,{boxcon, eq1, eq2, eq3});