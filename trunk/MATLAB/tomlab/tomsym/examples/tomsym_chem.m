%% Chemical Equilibrium Problem
% TomSym implementation of GAMS Example (CHEM,SEQ=21)
%
% The problem is to find the equilibrium composition of a
% mixture of different atoms.
%
% Bracken, J, and McCormick, G P, Chapter 5. In Selected Applications of
% Nonlinear Programming. John Wiley and Sons, New York, 1968, pp. 48-49.
%
% c: compounds (H, H2, H2O, N, N2, NH, NO, O, O2, OH)
%
% i: atoms (H  hydrogen, N  nitrogen, O  oxygen)

% Parameters

% Atoms per compound
a = [1   2    2  0   0   1   0  0   0   1;
    0   0    0  1   2   1   1  0   0   0;
    0   0    1  0   0   0   1  1   2   1];

% Number of elements in mixture
h=2; n=1; o=1;
mix = [h n o]';

% Gibbs free energy at 3500 k and 750 psi
gibbs = [-10.021; -21.096; -37.986; -9.846; -28.653; ...
    -18.918; -28.032; -14.640; -30.594; -26.11];

% Gibbs energy plus pressure
gplusp = gibbs + log(750*0.07031);

% Variables

% Number of mols in mixture
toms 10x1 x

% Component definition
eq1 = {};
for i=1:length(mix)
    eq1 = {eq1; sum(a(i,:)*x) == mix(i)};
end

% Energy definition (xb = total number of mols in mixture)
xb = sum(x);

% Total free energy in mixture
energy = sum(x.*(gplusp + log((x)/xb)));

% Bounds on x and xb
cbnd = {0.001 <= x; 0.01 <= xb};

solution = ezsolve(energy,{eq1, cbnd});