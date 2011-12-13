%% Flight Connections at a Hub
%
%% Problem description
% The airline SafeFlight uses the airport Roissy-Charles-de-Gaulle
% as a hub to minimize the number of flight connections to European
% destinations. Six Fokker 100 airplanes of this airline from
% Bordeaux, Clermont-Ferrand, Marseille, Nantes, Nice, and Toulouse
% are landing between 11am and 12:30pm. These aircraft leave for
% Berlin, Bern, Brussels, London, Rome, and Vienna between 12:30 pm
% and 13:30 pm. The numbers of passengers transferring from the
% incoming flights to one of the outgoing flights are listed in the
% table below.
%
% Numbers of passengers transferring between the different flights
%
%  +----------------+---------------------------------------+
%  |Origins         |          Destinations                 |
%  |                +------+----+--------+------+----+------+
%  |                |Berlin|Bern|Brussels|London|Rome|Vienna|
%  +----------------+------+----+--------+------+----+------+
%  |Bordeaux        |  35  | 12 |   16   |   38 |  5 |   2  |
%  |Clermont-Ferrand|  25  |  8 |    9   |   24 |  6 |   8  |
%  |Marseille       |  12  |  8 |   11   |   27 |  3 |   2  |
%  |Nantes          |  38  | 15 |   14   |   30 |  2 |   9  |
%  |Nice            |   –  |  9 |    8   |   25 | 10 |   5  |
%  |Toulouse        |   –  |  – |    –   |   14 |  6 |   7  |
%  +----------------+------+----+--------+------+----+------+
%
% For example, if the flight incoming from Bordeaux continues on to
% Berlin, 35 passengers and their luggage may stay on board during
% the stop at Paris. The flight from Nice arrives too late to be
% re-used on the connection to Berlin, the same is true for the
% flight from Toulouse that cannot be used for the destinations
% Berlin, Bern and Brussels (the corresponding entries in the table
% are marked with '–').
%
% How should the arriving planes be re-used for the departing
% flights to minimize the number of passengers who have to change
% planes at Roissy?
%
%% Variables
%
%  origdest                   People remaining in plane at transfer
%  idximp                     Impossible combinations
%
%% Reference
%
% Applications of optimization... Gueret, Prins, Seveaux

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Oct 7, 2005.   Last modified Apr 8, 2009.

%% Problem setup
origdest = [ 35 12 16 38  5  2;...
    25  8  9 24  6  8;...
    12  8 11 27  3  2;...
    38 15 14 30  2  9;...
    0  9  8 25 10  5;...
    0  0  0 14  6  7];

idximp = find(origdest == 0);

n = size(origdest,1);

cont = tom('cont',n,n,'int');

% All variables are integers.
bnds = {0 <= cont <= 1, cont(idximp) == 0};

% Orig constr.
con1 = {sum(cont,1) == 1};

% Dest constr.
con2 = {sum(cont,2) == 1};

% Objective
objective = -sum(sum(origdest.*cont));

constraints = {bnds, con1, con2};
options = struct;
options.solver = 'cplex';
options.name   = 'Flight Connections at a Hub';
sol = ezsolve(objective,constraints,[],options);

PriLev = 1;
if PriLev > 0
    p         = 6;                            % the number of planes
    in_names  = ['Bord';'Cler';'Mars';...     % city names in
        'Nant'; 'Nice'; 'Toul'];
    out_names = ['Berl';'Bern';'Brus';...     % city names out
        'Lond';'Rome';'Vien'];
    transfers = sol.cont;

    for i = 1:p,      % in  city has number i
        for o = 1:p,    % out city has number o
            if transfers(i,o) == 1,
                disp(['plane from '     in_names(i,:)  ...
                    ' continues to ' out_names(o,:) ])
            end
        end
    end
end

% MODIFICATION LOG
%
% 051021 med   Created.
% 060113 per   Added documentation.
% 060125 per   Moved disp to end
% 060203 med   Print level updated
% 090308 med   Converted to tomSym