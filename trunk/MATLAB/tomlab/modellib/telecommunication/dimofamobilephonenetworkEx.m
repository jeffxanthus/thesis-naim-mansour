% function Result = dimofamobilephonenetworkEx(PriLev)
%
% Creates a TOMLAB MIP problem for dimensioning of a mobile phone network
%
% DIMENSIONING OF A MOBILE PHONE NETWORK
%
% The figure below represents the typical architecture of a mobile
% phone network. Every elementary geographical zone or cell is served
% by a transmitter-receiver called a relay. The calls originating
% from a mobile phone first pass through these relays. Every relay is
% connected by cable or electro-magnetic waves to a transit node
% (hub). One of the hubs controls the network, this is the MTSO
% (Mobile Telephone Switching Office). A very reliable ring of fiber
% optic cable connects the hubs and the MTSO with high capacity
% links. It is capable of re-establishing itself in the case of a
% breakdown (self-healing ring) and there is no need to replicate it.
%
% At the present state of technology, there are no dynamic
% connections between the relays and the MTSO. The connections are
% fixed during the design phase, so it is necessary to choose the
% nodes of the ring that a relay should be connected to. The number
% of links between a cell c and the ring is called the diversity of
% the cell c, denoted by CNCTc. A diversity larger than 1 is
% recommended for making the system more reliable.
%
% The traffic in this kind of system is entirely digitized, expressed
% in values equivalent to bidirectional circuits at 64kbps (kilobit
% per second). This capacity corresponds to the number of
% simultaneous calls during peak periods. The ring has edges of known
% capacity CAP. The traffic TRAFc from a cell c is divided into equal
% parts (TRAFc / CNCTc) among the connections to the ring. This
% traffic is transmitted via the ring to the MTSO, that then routes
% it to another cell or to a hub that serves as the interface to the
% fixed-line telephone network. A relay may be connected directly to
% the MTSO that also has the functions of an ordinary hub.
%
% We consider a network of 10 cells and a ring of 5 nodes with a
% capacity of CAP = 48 circuits. The MTSO is at node 5. The following
% table indicates the traffic, the required number of connections and
% the cost per connection in thousand $ per cell. For example, cell 1
% is connectable with node 1 for a cost of $15,000. Its diversity is
% 2, which means it must be connected to at least two nodes of the
% ring. Its traffic capacity is of 22 simultaneous circuits. The
% objective is to define the connections of the cells to the ring
% that minimize the connection costs whilst remaining within the
% capacity limits and satisfying the constraints on the number of
% connections.
%
% Structure of a mobile phone network
% 
%    Cell 1
%    2 connections
% 
%    Relay ---------\
%                    \
%     |               \
%     |                \
%     |                 |
%     V                 V
%                
%    hub2 ============ hub3 <-- Relay  
%                               cell 2
%     ||                ||      1 connection
%     ||                ||
%     ||                ||
%
%    hub1 === MTSO === hub4
%
% Connection costs, traffic and number of connections per cell
%
% +------------+--+--+--+--+--+--+--+--+--+--+
% |Cell        | 1| 2| 3| 4| 5| 6| 7| 8| 9|10|
% +------------+--+--+--+--+--+--+--+--+--+--+
% |Hub 1       |15| 9|12|17| 8| 7|19|20|21|25|
% |Hub 2       | 8|11| 6| 5|22|25|25| 9|22|24|
% |Hub 3       | 7| 8| 7| 9|21|15|21|15|14|13|
% |Hub 4       |11| 5|15|18|19| 9|20|18|16| 4|
% |Hub 5 (MTSO)|10|14|15|24| 6|17|22|25|20|11|
% +------------+--+--+--+--+--+--+--+--+--+--+
% |Traffic     |22|12|20|12|15|25|15|14| 8|22|
% +------------+--+--+--+--+--+--+--+--+--+--+
% |Connections | 2| 2| 2| 2| 3| 1| 3| 2| 2| 2|
% +------------+--+--+--+--+--+--+--+--+--+--+
%
% VARIABLES
%
% hub_mat                    Matrix describing the hubs
% traffic                    Traffic from cells
% connections                Possible connections per cell
% capacity                   Capacity
%
% RESULTS
%
% For an interpretation of the results, try:
% Result   = dimofamobilephonenetworkEx(2);
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
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Nov 8, 2005.   Last modified Nov 8, 2005.

function Result = dimofamobilephonenetworkEx(PriLev)

if nargin < 1
   PriLev = 1;
end

hub_mat     = [15  9 12 17  8  7 19 20 21 25;...
      8 11  6  5 22 25 25  9 22 24;...
      7  8  7  9 21 15 21 15 14 13;...
      11  5 15 18 19  9 20 18 16  4;...
      10 14 15 24  6 17 22 25 20 11]*1000;

traffic     = [22 12 20 12 15 25 15 14  8 22]';
connections = [ 2  2  2  2  3  1  3  2  2  2]';
capacity    = 48;

Prob = dimofamobilephonenetwork(hub_mat, traffic, connections,...
   capacity);
Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   [h,c] = size(hub_mat) ;
   temp  = reshape(Result.x_k,c,h);
   for cell = 1:c,
      idx    = find(temp(cell,:));
      disp(['cell ' num2str(cell) ' connects to hub(s) ' num2str(idx)])
   end
end

% MODIFICATION LOG
%
% 051108 med   Created.
% 060116 per   Added documentation.
% 060126 per   Moved disp to end
