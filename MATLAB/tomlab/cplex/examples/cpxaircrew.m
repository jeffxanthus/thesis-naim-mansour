% cpxaircrew.m
%
% Air-crew schedule generation, ToD = Tours of Duty
%
% Based on D.M.Ryan, Airline Industry, Encyclopedia of Operations Research
% and Management Science
%
% Problem formulated at Optimization Theory, Linkoping University, called run5a
%
% Subfunctions used (at the end of this file):
%
%    function [M,T,S,R] = generateToDs()
%    function DATA = sectordata()
%
% function cpxaircrew(DEFPARAM,CUT)
%
% DEFPARAM  If =1 use default parameters, presolve, cuts, Dual simplex, 
%           If =0 no presolve (default)
% CUT       Cut strategy. 0 = no cuts, 1 = cuts, 2 = aggressive cuts
%           -1 = default CPLEX strategy. Default -1.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2001-2007 by Tomlab Optimization Inc., $Release: 10.0.0$
% Written May 1, 2001.  Last modified Feb 21, 2007.

function cpxaircrew(DEFPARAM,CUT)

if nargin < 2
   CUT = [];
   if nargin < 1
      DEFPARAM = [];
   end
end

if isempty(DEFPARAM), DEFPARAM = 0; end
if isempty(CUT),      CUT = -1; end

%clear all
global RELAX

% M is the MODIFIED adjacent matrix
% T is the salary time for ONE person
% S is the stayover cost for ONE person
% R saves info about the ToDs

RELAX = 0;
[M,T,S,R] = generateToDs;

RELAX = 1;
[M_r T_r S_r R_r] = generateToDs;

missed = find(sum(M')'==0);

M_exp = M; T_exp = T; S_exp = S; R_exp = R; I_exp = [];
for i=1:length(missed)
   flight = missed(i);
   fidx = find(M_r(flight,:)>0);
   M_exp = [M_exp M_r(:,fidx)];
   T_exp = [T_exp T_r(fidx)];
   S_exp = [S_exp S_r(fidx)];
   R_exp = [R_exp R_r(:,fidx)];
   I_exp = [I_exp flight*ones(1,length(fidx))];
end

A = sign(M);
for i=1:length(I_exp)
   col = zeros(34,1);
   col(I_exp(i)) = 1;
   A = [A col];         
end   


% ********************************************************************
% ********************************************************************

m = size(A,1); % number of flights (constraints)
n = size(A,2); % problem dimension

c   = 2*T_exp*150000/(365*24) + 2*S_exp;
c   = c(:);
x_L = zeros(n,1);
x_U = ones(n,1);

b_L = ones(m,1);
b_U = 100*ones(m,1);

MIP = 1

if MIP
   % Setting IntVars as 1:n (=length(x_L)) gives a pure IP solution
   IntVars = [1:n];
else
   % Setting IntVars as empty gives LP solution
   IntVars = [];
end

IntVars

callback = []; PriLev = []; Prob = [];
PI        = []; SC       = []; SI     = []; sos1   = []; sos2 = [];

Prob.P=1;

if CUT == 0, CUT = -1; end

if DEFPARAM
   % Empty cpxControl (default values) will give fastest execution
   cpxControl = [];
else
   % Setting Cplex Parameters (will slow down the computations)
   % Change to Primal simplex for root and subsolutions
   % 0 default 1 = Primal 2 = Dual (default) 3 = Network 4 = Barrier
   % cpxControl.NODEALG    = 0;
   % cpxControl.SUBALG     = 0;
   % No PreSolve
   cpxControl.PREIND       = 0;
end

   % Change cut strategy
   cpxControl.CLIQUES    = CUT;
   cpxControl.COVERS     = CUT;
   cpxControl.DISJCUTS   = CUT;
   cpxControl.FLOWCOVERS = CUT;
   cpxControl.FLOWPATHS  = CUT;
   cpxControl.FRACCUTS   = CUT;
   cpxControl.GUBCOVERS  = CUT;
   cpxControl.IMPLED     = CUT;
   cpxControl.MIRCUTS    = CUT;

cpxControl

% Cut strategy 0 (no cuts), 1 or 2 (Aggressive)

PriLev=1;

if 1

disp('Run problem on inequality form')
% INEQUALITY FORM
tic
[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes] = ...
    cplex(c, A, x_L, x_U, b_L, b_U, cpxControl, callback, PriLev, Prob, ...
           IntVars, PI, SC, SI, sos1, sos2);
toc
        
ix=find(x(1:n));
xprinti(ix,'x > 0');
xprint(x(ix),'xOpt');
xprinti(slack,'slack');

cpxPrint(PriLev,Inform,x,f_k,glnodes,lpiter,Prob,slack,v,basis,rc)

disp('End inequality form')
disp('Press return to continue')
pause
end

disp('Run problem on standard form')
% STANDARD FORM
if 1
   A = [A -eye(m)];
   c = [c;zeros(m,1)];
   b_U = b_L;
   %x_L = [x_L;zeros(m,1)];
   x_L = zeros(n+m,1);
   %x_U = [x_U;1000000*ones(m,1)];
   x_U = [ones(n,1);1000000*ones(m,1)];
else
   A = [A -eye(m) zeros(m,m);A zeros(m,m) eye(m)];
   c = [c;zeros(2*m,1)];
   b_L = [b_L;b_U];
   b_U = b_L;

   x_L = [x_L;zeros(2*m,1)];
   x_U = [x_U;1000*ones(2*m,1)];
end
if MIP
   % Setting IntVars as n (=length(x_L)) gives a pure IP solution
   IntVars = length(x_L);
else
   % Setting IntVars as empty gives LP solution
   IntVars = [];
end

tic
[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes] = ...
    cplex(c, A, x_L, x_U, b_L, b_U, cpxControl, callback, PriLev, Prob, ...
           IntVars, PI, SC, SI, sos1, sos2);
toc
        
ix=find(x(1:n));
xprinti(ix,'x > 0');
xprint(x(ix),'xOpt');
xprinti(slack,'slack');

cpxPrint(PriLev,Inform,x,f_k,glnodes,lpiter,Prob,slack,v,basis,rc)

% ********************************************************************
% generateToDs.m generates feasible ToDs for the Air-crew
% scheduling problem.
% ********************************************************************
%

function [M,T,S,R] = generateToDs()

M = []; % adjacent matrix
T = []; % salary time
S = []; % stayover cost
R = []; % saves info about the ToDs

% flight departure/arrival data
global DATA
DATA = sectordata;

for flight = 1:8
   
   ToD = zeros(34,1);
   ToD(flight) = 1;
   
   % Total stayover cost during ToD
   staycost = 0;
   
   % ToDtime is the total ToD-time
   ToDtime  = DATA(flight,7)-DATA(flight,4) + ...
              24*daydiff(DATA(flight,6),DATA(flight,3));
   
   % DUTYtime is the current time on duty
   DUTYtime = ToDtime;
   
   % DUTYjobs is the current number of flights on duty 
   DUTYjobs = 1;
   
   route = zeros(15,1);
   route(1) = flight;
   
   [M2 T2 S2 R2] = findToD(ToD,staycost,ToDtime,DUTYtime,DUTYjobs,route,flight);
   
   M = [M M2]; T = [T T2]; S = [S S2]; R = [R R2];
   
end
fprintf('\n  %d feasible ToD:s found\n',size(M,2));


% --------------------------------------
%         Subfunction findToD
% --------------------------------------
function [Mtmp,Ttmp,Stmp,Rtmp] = findToD(ToD,staycost,ToDtime,DUTYtime,DUTYjobs,route,flight)

global DATA RELAX

STAYCOSTS = [219 250 322 180 181 315]/24; % stayovercosts
AKL  = 1;

if ToDtime > 5*24 % infeasible ToD
   Mtmp = []; Ttmp = []; Stmp = []; Rtmp = [];
   disp('  TO LONG TRIP !!');
   return
end

if DATA(flight,5)==AKL % if back in Auckland
   Mtmp = ToD;
   Ttmp = ToDtime;
   Stmp = staycost;
   Rtmp = route;
   return
end

current_city = DATA(flight,5); 
current_day  = DATA(flight,6);
current_time = DATA(flight,7);

% find all flights with current city as departure.
tmpidx = find(DATA(:,2)==current_city);

idx = [];
for i=1:length(tmpidx)
   flag = 1;
   flight_i = tmpidx(i);
   if ~( (daydiff(DATA(flight_i,3),current_day)==0) & ...
         (DATA(flight_i,4) < current_time) ) & ...
         (daydiff(DATA(flight_i,3),current_day)<6)
      if ToDtime + DATA(flight_i,7)-current_time +...
         24*daydiff(DATA(flight_i,6),current_day) > 24*5
         flag = 0;
      end
      if RELAX==0
         if DATA(flight_i,4)-current_time+24*daydiff(DATA(flight_i,3),current_day) < 16 & ...
            DUTYtime + DATA(flight_i,7)-current_time + 24*daydiff(DATA(flight_i,6),current_day) > 12
            flag = 0;
         end
         if DATA(flight_i,4)-current_time+24*daydiff(DATA(flight_i,3),current_day) < 16 & ...
            DUTYjobs >= 2
            flag = 0;
         end
      end
      if flag
         idx = [idx flight_i];
      end
   end
end
% idx now contains all feasible flights but rest may be necessary


Mtmp = []; Ttmp = []; Stmp = []; Rtmp = [];
for i=1:length(idx)
   
   flight_i = idx(i);
   
   ToD_i = ToD;
   ToD_i(flight_i) = max(ToD)+1;
   
   % time from arrival to departure
   t1 = DATA(flight_i,4)-current_time + ... 
        24*daydiff(DATA(flight_i,3),current_day);
     
   % time of new flight  
   t2 = DATA(flight_i,7)-DATA(flight_i,4) + ...
        24*daydiff(DATA(flight_i,6),DATA(flight_i,3));
   
   ToDtime_i  = ToDtime + t1 + t2;
   
   route_i = route;
   
   if t1 >= 16 % must take a break
      DUTYtime_i = t2;
      DUTYjobs_i = 1;
      route_i(sum(sign(route_i))+1) = 1000;
   else
      DUTYtime_i = DUTYtime + t1 + t2;
      DUTYjobs_i = DUTYjobs + 1;
   end
   
   route_i(sum(sign(route_i))+1) = flight_i;
   
   staycost_i = staycost + t1*STAYCOSTS(current_city);
   
   [Mtmp_i Ttmp_i Stmp_i Rtmp_i] = ...
      findToD(ToD_i,staycost_i,ToDtime_i,DUTYtime_i,DUTYjobs_i,route_i,flight_i);
   
   Mtmp = [Mtmp Mtmp_i];
   Ttmp = [Ttmp Ttmp_i];
   Stmp = [Stmp Stmp_i];
   Rtmp = [Rtmp Rtmp_i];
      
end


% -----------------------------------

function days = daydiff(day_dep,day_arr)
if day_dep >= day_arr
   days = day_dep-day_arr;
else
   days = 7-day_arr+day_dep;
end

% ********************************************************************
% function sectordata() returns sector data with departure
% and arrival times in hours converted to GMT.
% ********************************************************************

function DATA = sectordata()

% Assign a number to each city
AKL = 1; LAX = 2; HNL = 3; APW = 4; TBU = 5; NAN = 6;

% Assign a number to each day of week
MON = 1; TUE = 2; WED = 3; THU = 4; FRI = 5; SAT = 6; SUN = 7;

DATA = [...
1	AKL	MON	2145	NAN	TUE	50
2	AKL	TUE	1330	TBU	TUE	1715
3	AKL	WED	1145	HNL	TUE	2230
4	AKL	THU	1145	HNL	WED	2230
5	AKL	FRI	1000	NAN	FRI	1305
6	AKL	SAT	1145	HNL	FRI	2230
7	AKL	SAT	1500	NAN	SAT	1805
8	AKL	SUN	1500	NAN	SUN	1805
9	LAX	MON	2130	HNL	TUE	30
10	LAX	SAT	2130	HNL	SUN	30
11	LAX	SUN	2300	HNL	MON	200
12	HNL	MON	1130	LAX	MON	2005
13	HNL	SAT	450	LAX	SAT	1325
14	HNL	SUN	450	LAX	SUN	1325
15	HNL	TUE	150	NAN	WED	640
16	HNL	TUE	445	APW	TUE	910
17	HNL	TUE	2345	AKL	THU	705
18	HNL	WED	2345	NAN	FRI	435
19	HNL	THU	2345	AKL	SAT	705
20	HNL	FRI	2345	AKL	SUN	705
21	HNL	SUN	150	NAN	MON	640
22	HNL	MON	320	NAN	TUE	810
23	APW	MON	2100	HNL	TUE	315
24	APW	TUE	1025	TBU	WED	1155
25	TBU	WED	1300	AKL	WED	1505
26	TBU	TUE	1815	APW	MON	1945
27	NAN	TUE	150	HNL	MON	1015
28	NAN	FRI	1405	HNL	THU	2230
29	NAN	SAT	1905	HNL	SAT	330
30	NAN	SUN	1905	HNL	SUN	330
31	NAN	WED	740	AKL	WED	1045
32	NAN	FRI	535	AKL	FRI	840
33	NAN	MON	740	AKL	MON	1045
34	NAN	TUE	910	AKL	TUE	1215 ];

% convert from hours and minutes to hours
for i=1:34
   tmp1 = floor(DATA(i,4)/100);
   DATA(i,4) = tmp1+(DATA(i,4)-100*tmp1)/60;
   tmp2 = floor(DATA(i,7)/100);
   DATA(i,7) = tmp2+(DATA(i,7)-100*tmp2)/60;
end


% convert to GMT
for i=1:34
   for j=2:3:5
      
      % LAX, HNL and APW
      if DATA(i,j)==LAX | DATA(i,j)==HNL | DATA(i,j)==APW
         if DATA(i,j) == LAX
            corr = 8;
         elseif DATA(i,j) == HNL
            corr = 10;
         else
            corr = 11;
         end
         if 24-DATA(i,j+2) > corr
            DATA(i,j+2) = DATA(i,j+2) + corr;
         else
            DATA(i,j+2) = corr - (24-DATA(i,j+2));
            if DATA(i,j+1) == 7 % if sunday
               DATA(i,j+1) = 1;
            else
               DATA(i,j+1) = DATA(i,j+1) + 1;
            end
         end
      end
      
      % AKL, TBU and NAN
      if DATA(i,j)==AKL | DATA(i,j)==TBU | DATA(i,j)==NAN
         if DATA(i,j) == NAN
            corr = 12;
         else
            corr = 13;
         end
         
         if DATA(i,j+2) > corr
            DATA(i,j+2) = DATA(i,j+2) - corr;
         else
            DATA(i,j+2) = 24-(corr-DATA(i,j+2));
            if DATA(i,j+1) == 1 % if sunday
               DATA(i,j+1) = 7;
            else
               DATA(i,j+1) = DATA(i,j+1) - 1;
            end
         end
      end
      
   end
end


if 1
   fprintf('\n   Flight ID        Departure                  Arrival');
   fprintf('\n     ID         City   Day    Time        City   Day    Time\n');
   for i=1:34
      fprintf('\n %6d       %5d %5d   %7.3f     %5d %5d   %7.3f',...
         i,DATA(i,2),DATA(i,3),DATA(i,4),DATA(i,5),DATA(i,6),DATA(i,7));
   end
   fprintf('\n\n\n');
end

% MODIFICATION LOG:
%
% 050113 med Modification log added
% 050117 med [] removed for x_0
% 050209 med Changed name from aircrew to cpxaircrew
% 070221 hkh Revise IntVars format
