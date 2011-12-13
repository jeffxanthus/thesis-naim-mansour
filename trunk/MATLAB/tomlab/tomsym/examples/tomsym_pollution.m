%% Industrial Pollution Control
% TomSym implementation of GAMS Example (POLLUT,SEQ=96)
%
% This is an example of planning production within certain pollution
% standards.
%
% Mangasarian, O L, Nonlinear Programming. McGraw Hill, New York, 1973.
%
% j: (food, textile, apparel, lumber, furniture, pulp+paper
%     printing, chemicals, coal+petro, rubber, leather
%     clay+stone, steel, nf-metals, metals-pr, machinery
%     elec-mach, transport, prec-mach, misc)
%
% i:  pollutant (cod, so2, land, water)

% Matrix w(j,i) unit load
w = [.08488  .00909  .02990  .07195;
    .03353  .02470  .08980  .17032;
    .03353  .02470  .02780  .01345;
    .00135  .00084  .06890  .01606;
    .00153  .00084  .03500  .01429;
    .23368  .07461  .04690  .15532;
    .07609  .05767  .01290  .01850;
    .07689  .05767  .05370  .09680;
    .03736  .01663  .27080  .02900;
    .02794  .00456  .06040  .09200;
    .02794  .00456  .04980  .07320;
    .00213  .08800  .11510  .09520;
    .00633  .02361  .08900  .07780;
    .00091  .03376  .02670  .03720;
    .00125  .00860  .06790  .04930;
    .00089  .00376  .06210  .02480;
    .00123  .00369  .02950  .02300;
    .00079  .00127  .07870  .04300;
    .00396  .00252  .03830  .03010;
    .00931  .00252  .04240  .03500];

% Matrix z(j,*) other data
z =  [9.600   .121   .1064    29406   24709;
    6.353   .194   .1611    21375   18918;
    9.818   .204   .0848     8423   20636;
    7.371   .181   .1428    13873   10457;
    10.220   .171   .0955     8470    9739;
    6.255   .145   .1759    36375   18880;
    8.149   .304   .2190    66016   44480;
    7.794   .146   .1439    80134   36526;
    8.400   .173   .1314     1327     758;
    9.933   .174   .0987     4414    4921;
    11.069   .167   .0817     3709    6766;
    6.528   .192   .1665    14496    9368;
    7.928   .116   .1448   102399   31127;
    10.559   .091   .1026    28008    1166;
    6.606   .227   .1567    69314   59525;
    7.153   .208   .1506    90014   63048;
    11.146   .151   .0895    29360   29839;
    6.884   .199   .1602    27687   16945;
    6.660   .253   .1446     4009    4828;
    7.929   .182   .1256    24323   24569];

% Scalars
a = .6;
b = 1.4;
gamma1 = 1.4;
gamma2 = 0.9;

% Parameter tau(i)
tau = [.153e6 .12e6 .25e6 .25e6];

% Variables (k = capital, l = labor)
k = tom('k',size(w,1),1);
l = tom('l',size(w,1),1);

output =  sum(z(:,1).*k.^(1-z(:,2)).*l.^z(:,2));

tk = sum(k);
tl = sum(l);

constr = {
    sum(w.*repmat(k./z(:,3),1,size(w,2))) <= tau
    tk >= gamma2*tl
    tk <= gamma1*tl
    a*z(:,4) <= k <= b*z(:,4)
    a*z(:,5) <= l <= b*z(:,5)
    };

guess = [];

options = struct;
options.name = 'pollution';
solution = ezsolve(-output,constr,guess,options);
