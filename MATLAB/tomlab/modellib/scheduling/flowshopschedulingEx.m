% function Result = flowshopschedulingEx(PriLev)
%
% Creates a TOMLAB MIP problem for flow shop scheduling
%
% FLOW-SHOP SCHEDULING
%
% A workshop that produces metal pipes on demand for the automobile
% industry has three machines for bending the pipes, soldering the
% fastenings, and assembling the links. The workshop has to produce
% six pieces, for which the durations of the processing steps (in
% minutes) are given in the following table. Every workpiece first
% goes to bending, then to soldering, and finally to assembly of the
% links. Once started, any operations must be carried out without
% interruption, but the workpieces may wait between the machines.
%
% Processing durations in minutes
%
% +---------+-+-+-+-+-+-+
% |Workpiece|1|2|3|4|5|6|
% +---------+-+-+-+-+-+-+
% |Bending  |3|6|3|5|5|7|
% |Soldering|5|4|2|4|4|5|
% |Assembly |5|2|4|6|3|6|
% +---------+-+-+-+-+-+-+
%
% Every machine only processes one piece at a time. A workpiece may
% not overtake any other by passing onto the following machine. This
% means that if at the beginning a sequence of the workpieces is
% established, they will be processed on every machine in exactly
% this order. Which is the sequence of workpieces that minimizes the
% total time for completing all pieces?
%
% VARIABLES
%
% nummachines                The number of machines
% proctimes                  Time to process a workpiece in a machine
%
% RESULTS
%
% The vector x_k can be interpreted by running the following:
% Result = flowshopschedulingEx(2);
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
% Written Oct 13, 2005.   Last modified Oct 13, 2005.

function Result = flowshopschedulingEx(PriLev)

if nargin < 1
   PriLev = 1;
end

nummachines = 3;
proctimes   = [3 6 3 5 5 7;...
      5 4 2 4 4 5;...
      5 2 4 6 3 6];

Prob = flowshopscheduling(nummachines, proctimes);

Result = tomRun('cplex', Prob, PriLev);

if PriLev > 1,
   p      = size(proctimes,2);                % number of pieces
   m      = nummachines;                      % number of machines
   order  = reshape(Result.x_k(1:p*p),p,p);   % the order of the pieces
   order(find(order<0.5)) = 0;                % delete bad zeros
   wait1  = Result.x_k(p*p+1:p*p+p);          % wait before machine 1
   wait2  = Result.x_k(end-2*p+1:end-p);      % wait before machine 2
   wait3  = Result.x_k(end-p+1:end);          % wait before machine 3
   proctimes   = [3 6 3 5 5 7;...             % the processing times
         5 4 2 4 4 5;...
         5 2 4 6 3 6];
   sequence = [] ;                            % the sequence of pieces
   for i = 1:p,
      sequence = [sequence find(order(:,i))]; % insert piece by piece
   end
   disp(['A best sequence is: ' num2str(sequence)])
   
   mt1 = [];   % empty machine-times for machine 1
   mt2 = [];   % empty machine-times for machine 2
   mt3 = [];   % empty machine-times for machine 3
   t11 = 0;    % timepoint for entering machine 1
   t12 = 0;    % timepoint for exiting  machine 1
   t21 = 0;    % timepoint for entering machine 2
   t22 = 0;    % timepoint for exiting  machine 2
   t31 = 0;    % timepoint for entering machine 3
   t32 = 0;    % timepoint for exiting  machine 3
   first2 = 1; % the first piece in machine 2
   first3 = 1; % the first piece in machine 3
   
   for i = 1:length(sequence),
      id  = sequence(i);
      
      % times in and out of machine 1
      t11  = t12 + wait1(id);       
      t12  = t11 + proctimes(1,id); 
      
      % times in and out of machine 2
      if first2 == 1,
         t21 = t12 + wait2(id);
         first2 = 0 ;
      else
         t21 = t22 + wait2(id);
      end
      t22  = t21 + proctimes(2,id);
      
      % times in and out of machine 3
      if first3 == 1,
         t31 = t22 + wait3(id);
         first3 = 0 ;
      else
         t31 = t32 + wait3(id);
      end
      t32  = t31 + proctimes(3,id);
      
      % times in and out of machines
      mt1 = [mt1; t11  t12];
      mt2 = [mt2; t21  t22];
      mt3 = [mt3; t31  t32];
   end
   
   mt = [mt1 mt2 mt3];
   
   disp('FLOW FOR THE PIECES')
   for j = 1:p,
      disp(['piece ' num2str(sequence(j)) ' has this flow' ])
      for k = 1:m,
         disp(['   machine ' num2str(k) ': ' num2str(mt(j,(k-1)*2+1)) '-' num2str(mt(j,(k-1)*2+2)) ])
      end
   end
   disp('FLOW FOR THE MACHINES')
   for k = 1:m,
      disp(['machine ' num2str(k) ' has this flow' ])
      for l = 1:p,
         disp(['   piece ' num2str(sequence(l)) ': ' ... 
               num2str(mt(l,((k-1)*2+1))) '-' num2str(mt(l,((k-1)*2+2)))])
      end
   end
   
end

% MODIFICATION LOG
%
% 051010 med   Created.
% 060111 per   Added documentation. 
% 060124 per   Interpretation of results upgraded.
% 060126 per   Moved disp to end
