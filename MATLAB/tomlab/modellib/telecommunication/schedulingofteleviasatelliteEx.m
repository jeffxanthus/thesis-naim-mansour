% function Result = schedulingofteleviasatelliteEx(PriLev)
%
% Creates a TOMLAB MIP problem for scheduling of telecommunications via 
% satellite
%
% SCHEDULING OF TELECOMMUNICATIONS VIA SATELLITE
%
% A digital telecommunications system via satellite consists of a
% satellite and a set of stations on earth which serve as interfaces
% with the terrestrial network. With SS-TDMA (Satellite-Switch, Time
% Division Multiple Access) access technology, the satellite divides
% its time among the stations. Consider for example the transmissions
% from four transmitter stations in the US to four receiver stations
% in Europe. The following table gives a possible 4 × 4 traffic
% matrix. TRAFtr is the quantity of data transmitted from station t
% to station r. It is expressed in seconds of transmission duration,
% because all lines have the same constant transmission rate.
%
% Matrix TRAF with its lower bound LB
%
%  +----+--+--+--+--+-------+ 
%  |TRAF| 1| 2| 3| 4| rowt  |
%  +----+--+--+--+--+-------+ 
%  |  1 | 0| 7|11|15|   33  |
%  |  2 |15| 8|13| 9|   45  |
%  |  3 |17|12| 6|10|   45  |
%  |  4 | 6|13|15| 4|   38  |
%  +----+--+--+--+--+-------+ 
%  |colr|38|40|45|38|LB = 45|
%  +----+--+--+--+--+-------+ 
%
% The satellite has a switch that allows any permutation between the
% four transmitters and the four receivers. The figure below gives an
% example of a permutation connecting the transmitters 1 to 4 to the
% receivers 3, 4, 1, and 2 respectively. These connections allow
% routing a part of the traffic matrix, called a mode. A part of a
% matrix entry transmitted during a mode is a packet of data. A mode
% is thus a 4 × 4 matrix M with at most one non-zero packet per row
% and column and such that Mtr <= TRAFtr for all t, r. To every
% mode corresponds a slice of a schedule as shown in the figure.
%
%  +----+--+--+--+--+         +---------+---------------+
%  |TRAF| 1| 2| 3| 4|         |Stations | Packets       |
%  +----+--+--+--+--+         +---------+---------------+
%  |  1 | 0| 0|11| 0|         | 1 --> 3 |-----11--->    |
%  |  2 | 0| 0| 0| 9|         | 2 --> 4 |------9->      |
%  |  3 |15| 0| 0| 0|         | 3 --> 1 |-----15------->|
%  |  4 | 0|13| 0| 0|         | 4 --> 2 |-----13----->  |
%  +----+--+--+--+--+         +---------+---------------+
%  |colr|38|40|45|38|LB = 45
%  +----+--+--+--+--+ 
%
%
% Example of a mode and associated schedule
%  
% A valid schedule of transmissions defines a sequence of
% permutations for the on-board switch that routes all the traffic of
% the matrix TRAF. In other words, this boils down to decomposing
% TRAF into a sum of mode matrices. An element of TRAF may be
% fragmented, like TRAF31 that is only partially transmitted in the
% mode represented in the figure above. A fragmented element will
% appear in several packets and several modes of the final
% decomposition. The duration of a mode is the length of its longest
% packet. The objective is to find a schedule with minimal total duration.
%
% VARIABLES
%
% traffic                    A matrix describing the traffic
%
% RESULTS
%
% For an interpretation of the results, let PriLev > 1,
% schedulingofteleviasatelliteEx(2);
%
% REFERENCES
%
% Applications of optimization... Gueret, Prins, Seveaux
% http://web.univ-ubs.fr/lester/~sevaux/pl/index.html
%
% INPUT PARAMETERS
% PriLev       Print Level
%
% OUTPUT PARAMETERS
% Result       Result structure.

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2007 by Tomlab Optimization Inc., $Release: 6.0.0$
% Written Nov 22, 2005.   Last modified Dec 18, 2007.

function Result = schedulingofteleviasatelliteEx(PriLev)

if nargin < 1
   PriLev = 1;
end

if PriLev > 1,    % if:   PriLev > 1
   sequence = []; % then: let us catch all x_k 
end

traffic      = [  0   7  11  15;...
      15   8  13   9;...
      17  12   6  10;...
      6  13  15   4];

TQBS = schedulingofteleviasatellite(traffic);

while sum(sum(TQBS)) > 0.01
   n1   = size(TQBS,1);
   n2   = n1^2;
   n    = n2 + 1;
   x_L  = zeros(n,1);
   x_U  = [ones(n2,1);inf];
   IntVars = [ones(n2,1);0];
   b_L1 = ones(n1,1);
   b_U1 = ones(n1,1);
   b_L2 = ones(n1,1);
   b_U2 = ones(n1,1);
   b_L3 = zeros(n1,1);
   b_U3 = inf*ones(n1,1);
   A1   = zeros(n1,n);
   A2   = zeros(n1,n);
   A3   = zeros(n1,n);
   for i=1:n1
      idx1  = [(i-1)*n1+1:i*n1];
      idx2  = find(TQBS(:,i) == 0);
      idx1(idx2) = [];
      idx3  = [i:n1:n2-n1+i];
      idx4  = find(TQBS(i,:) == 0);
      idx3(idx4) = [];
      A1(i,idx1) = ones(1,length(idx1));
      A2(i,idx3) = ones(1,length(idx3));
      A3(i,(i-1)*n1+1:i*n1) = TQBS(:,i)';
      A3(i,end) = -1;
   end
   c = [zeros(n2,1);-1];
   Prob = mipAssign(c, [A1;A2;A3], [b_L1;b_L2;b_L3], [b_U1;b_U2;b_U3], x_L, x_U, [], 'Satellite Scheduling',...
      [], [], IntVars);
   Result = tomRun('cplex', Prob, 0);
   x_k = Result.x_k;
   x_k = reshape(x_k(1:end-1,1), n1, n1); 
   TQBS = TQBS - x_k*(-Result.f_k);
   
   if PriLev > 1,
      sequence = [sequence Result.x_k];
   end
   
end
PrintResult(Result,PriLev);

if PriLev > 1,
   for t = 1:size(sequence,2),
      disp(['Transmission ' num2str(t) ' transfers ' num2str(sequence(end,t)) ' packet(s) of data' ])
      this = reshape(sequence(1:end-1,t),n1,n1);
      this(find(this < 0.5))  = 0; % remove bad zeros
      for j = 1:n1,
         disp(['   from station ' num2str(j) ' to ' num2str(find(this(j,:))) ])
      end
   end
end

% MODIFICATION LOG
%
% 051122 med   Created.
% 060116 per   Added documentation.
% 060126 per   Moved disp to end
% 060131 per   cleaned up help text
% 071218 ango  Multiple CPLEX solutions handled gracefully