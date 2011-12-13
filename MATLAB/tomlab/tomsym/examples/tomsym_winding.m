%% Winding Factor of Electrical Machines
% TomSym implementation of GAMS Example (WINDFAC, SEQ=224)
%
% This model determines the optimal winding factor for electrical
% machines.
%
% Michna, M, and Gdanska, P, Winding Factor of Electrical Machines, 1984.

clear mex % minlpBB sometimes needs this

ms = 3; % number of phases
p  = 2; % number of pole pairs
K  = 5; % harmonic order
ns = 1; % coil span
clear pi

%     q:  number of slots per one phase and per one pole
%     Nz: number of slots
%     s:  span

toms integer q

% Variable bounds
cbnd = {1 <= q <= 10};

toms kz1 kz3 kz5 ks1 ks3 ks5 % coil-group factor (1-3) coil-span factor
toms kw kw3 kw5 % winding factor
toms kw1 % inding factor for first harmonic

% kw1 bound
cbnd = {cbnd; 0.8 <= kw1};

Nz = 2*ms*q*p;

alfae = (2*pi*p)/Nz;

% Slots pitch
tauz = Nz/(2*p); % == ms*q

s = tauz - ns;

eqs = {
    (q * sin( alfae / 2)) * kz1 == sin(q * alfae / 2)
    ks1 == sin((s * pi) / (tauz * 2))
    kw1 == ks1 * kz1
    (q * sin(3 * alfae / 2)) * kz3 == sin(3 * q * alfae / 2)
    ks3 == sin((3 * s * pi) / (tauz * 2))
    kw3 == ks3 * kz3
    (q * sin(5 * alfae / 2)) * kz5 == sin(5 * q * alfae / 2)
    ks5 == sin((5 * s * pi) / (tauz * 2))
    kw5 == ks5 * kz5
    kw == kw3^2 + kw5^2
};

options = struct;
options.name = 'WindFac';
options.solver = 'KNITRO';
guess = {ks1 == 0.5, ks3 == 0.5, ks5 == 0.5};

solution = ezsolve(kw,{cbnd;eqs},guess,options);
