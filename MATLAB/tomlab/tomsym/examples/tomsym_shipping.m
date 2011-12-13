%% Structural Optimization
% TomSym implementation of GAMS Example (SHIP,SEQ=22)
%
% This model designs a vertically corrugated transverse bulkhead of an
% oil tanker. The objective is to design for minimum weight and
% meet stress, moment of inertia and plate thickness constraints.
%
% Bracken, J, and McCormick, G P, Chapter 6. In Selected Applications of
% Nonlinear Programming. John Wiley and Sons, New York, 1968.

% s: bulkhead sections (top, middle, bottom)

toms 3x1 s

gam = 0.001; % specific gravity of water
sig = 1200;  % maximum bending stress
dnv = 3.9;   % det norske veritas factor
ha = 250;    % height above panel
gamsteel = .0078; % specific weight of steel

% Length of panel
l = [495; 385; 315];

% Height at the base of panel
hb = ha+sum(full(spdiags(repmat(l',3,1),[0 1 2],3,3)))';

% Height at the middle of panel
h = hb - l/2;

% Constant number one
k1 = gam*h.*l.*l/12/sig;

% Constant number two
k2 = dnv*1.05e-4*sqrt(hb);

e = 0.8;
width = 500;
ca = 0.2;
tlow = 1.05;

toms 3x1 z th % module, plate thickness
% width of flange, length of web, depth of corrugation, width of
% corrugation, weight of structure
toms wl lw d wc w

% Module definition
eq1 = {z == d*th*(lw/3+wl*e)/2};

% Width of corrugation - definition
eq2 = {(wc-wl)^2 == lw^2 - d^2};

% Bending stress
eq3 = {z >= k1*wc};

% Moment of inertia
eq4 = {z*d/2 >= 2.2*(k1*wc).^(4/3)};

% Plate thickness - width of flange
eq5 = {th >= k2*wl + ca};

% Plate thickness - length of web
eq6 = {th >= k2*lw + ca};

% Geometric constraint
eq7 = {lw >= d};

% Total weight of structure

obj = gamsteel*width*(wl+lw)*sum(th.*l)/wc/1000;

% Bound on th
cbnd = {th >= tlow};

% Starting point
x0 = {th == [1.2;1.2;1.3]; wl == 45.8; lw == 43.2
    d == 30.5; wc == wl+sqrt(lw^2-d^2)
    z == d*th*(lw/3+wl*e)/2};

cons = {cbnd; eq1; eq2; eq3; eq4; eq5; eq6; eq7};
solution = ezsolve(obj,cons,x0);