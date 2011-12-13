%% Efficiency of Hospitals
%
%% Problem description
% The administration of the hospitals in Paris decides to measure the
% efficiency of the surgery departments in four major hospitals with
% a desire to improve the service to the public. To keep this study
% anonymous, the hospitals are named H1 to H4. The method suggested
% to measure the efficiency is DEA (Data Envelopment Analysis). This
% method compares the performance of a fictitious hospital with the
% performances of the four hospitals.
%
% Three initial indicators (resources) are taken into account: the
% number of non-medical personnel, the general expenses, and the
% available number of beds. In addition, four final indicators
% (services) are analyzed: the number of hospital admissions per day,
% the number of consultations in the outpatients clinic, the number
% of nurses on duty every day, and the number of interns and doctors
% on duty every day. The corresponding data have been analyzed over a
% period of two years and the numbers representing a day of average
% activity in every hospital are given in the following two tables.
%
% Resource indicators
%
%  +---------------------+-----+------+-----+-----+
%  |                     | H1  | H2   | H3  | H4  |
%  +---------------------+-----+------+-----+-----+
%  |Non-medical personnel| 90  | 87   | 51  | 66  |
%  |General expenses (k$)|38.89|109.48|40.43|48.41|
%  |Number of beds       | 34  | 33   | 20  | 33  |
%  +---------------------+-----+------+-----+-----+
%
% Service indicators
%
%  +-------------------+-----+-----+-----+-----+
%  |                   | H1  | H2  | H3  | H4  |
%  +-------------------+-----+-----+-----+-----+
%  |Admissions         |30.12|18.54|20.88|10.42|
%  |Consultations      |13.54|14.45| 8.52|17.74|
%  |Interns and doctors| 13  |  7  |   8 | 26  |
%  |Nurses on duty     | 79  | 55  |  47 | 50  |
%  +-------------------+-----+-----+-----+-----+
%
% Justify through the DEA method how hospital H2 is performing
% compared to the others.
%
%% Variables
%
%  resources                  A matrix describing the resources
%  services                   A matrix describing the services
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
resources  = [90 87 51 66;...
    38.89 109.48 40.43 48.41;...
    34 33 20 33];
services   = [30.12 18.54 20.88 10.42;...
    13.54 14.45 8.52 17.74;...
    13 7 8 26;...
    79 55 47 50];

h = size(resources,2);    %coef
s = size(services,1);     %services
r = size(resources,1);    %resources

coef = tom('coef',h,1);
fserv = tom('fserv',s,1);
fres = tom('fres',r,1);
eff = tom('eff',1,1);

% No variables are binary
bnds = {coef >= 0, fserv >= 0, fres >= 0, eff >= 0};

% Coef constraints
con1 = {sum(coef) == 1};

% Service constraint
con2 = {services*coef == fserv};

% Resource constraint
con3 = {resources*coef == fres};

% Objective
objective = eff;

indices = zeros(size(resources,2),1);
for i=1:size(resources,2)
    % Service indicators, greater than ficticious ones
    con4 = {fserv >= services(:,i)};

    % Efficiency relationship
    con5 = {fres <= resources(:,i)*eff};

    constraints = {bnds, con1, con2, con3, con4, con5};
    options = struct;
    options.solver = 'cplex';
    options.name   = 'Efficiency of Hospitals';
    sol = ezsolve(objective,constraints,[],options);

    indices(i,1) = sol.eff;
end

PriLev = 1;
if PriLev > 0
    for i = 1:length(indices),
        disp(['H' num2str(i) ' is ' num2str(100*indices(i)) '% efficient'])
    end
end

% MODIFICATION LOG
%
% 051206 med   Created
% 060118 per   Added documentation
% 060125 per   Moved disp to end
% 090325 med   Converted to tomSym