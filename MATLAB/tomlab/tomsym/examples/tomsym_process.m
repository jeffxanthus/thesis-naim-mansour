%% Alkylation Process Optimization
% TomSym implementation of GAMS Example (PROCESS,SEQ=20)
%
% Optimization of a alkylation process.
%
% Bracken, J, and McCormick, G P, Chapter 4. In Selected Applications
% of Nonlinear Programming. John Wiley and Sons, New York, 1968.
%
% olefin:   olefin feed (bpd)
%
% isor:     isobutane recycle (bpd)
%
% acid:     acid addition rate (1000lb per day)
%
% alkylate: alkylate yield (bpd)
%
% isom:     isobutane makeup (bpd)
%
% strength: acid strength (weight pct)
%
% octane:   motor octane number
%
% ratio:    isobutane makeup to olefin ratio
%
% dilute:   acid dilution factor
%
% f4:       f-4 performance number

toms olefin isor acid alkylate isom strength octane ratio dilute f4
x = [olefin isor acid alkylate isom strength octane ratio dilute f4]';

% Variables are positive
cbnd = {x >= 0};

% Alkylate yield definition
eq1 = {alkylate == olefin*(1.12+.13167*ratio-.00667*(ratio^2))};

% Isobutane makeup definition
eq2 = {1.22*alkylate == olefin + isom};

% Acid strength definition
eq3 = {acid == alkylate*dilute*strength/(98-strength)/1000};

% Motor octane number
eq4 = {octane == 86.35+1.098*ratio-.038*(ratio^2)-.325*(89-strength)};

% Isobutane to olefin ratio
eq5 = {ratio == (isor+isom)/olefin};

% Dilution definition
eq6 = {dilute == 35.82 - .222*f4};

% f-4 definition
eq7 = {f4 == -133 + 3*octane};

% Profit definition
profit = .063*alkylate*octane - 5.04*olefin -...
    .035*isor - 10*acid - 3.36*isom;

% Ranged alkylate yield definition
rangey = olefin*(1.12+.13167*ratio-.00667*(ratio^2))/alkylate;

% Ranged motor octane number
rangem = (86.35+1.098*ratio-.038*ratio^2-.325*(89-strength))/octane;

% Ranged dilution definition
ranged = (35.82 - .222*f4)/dilute;

% Ranged f-4 definition
rangef = (-133 + 3*octane)/f4;

% Bounds
eq8 = {0.9 <= rangey <= 1.1
    0.9 <= rangem <= 1.1
    0.9 <= ranged <= 1.1
    0.9 <= rangef <= 1.1};

eq9 = {85 <= strength <= 93
    90 <= octane <= 95
    3 <= ratio <= 12
    1.2 <= dilute <= 4
    145 <= f4 <= 162
    10 <= olefin <= 2000
    isor <= 16000
    acid <= 120
    alkylate <= 5000
    isom <= 2000};

% Starting point
x0 = {olefin == 1745; isor == 12000
    acid == 110; alkylate == 3048; isom == 1974;
    strength == 89.2; octane == 92.8;
    ratio == 8; dilute == 3.6; f4 == 145};

con1 = {cbnd; eq1; eq2; eq3; eq4; eq5; eq6; eq7; eq9};
solution1 = ezsolve(-profit,con1,x0);

con2 = {cbnd; eq2; eq3; eq5; eq8; eq9};
solution2 = ezsolve(-profit,con2,x0);