%
% cpxBuildConflict.m
%
% cpxBuildConflict provides a convenient method for generating conflict
% refinement groups for use with the Conflict Refinement feature of CPLEX
%
% SYNTAX:
%
% (1) function confgrps = cpxBuildConflict(Prob,mode)
% OR
% (2) function confgrps = cpxBuildConflict(n,m_lin,m_quad,m_sos,m_log,'mode')
%
% Inputs for (1): function confgrps = cpxBuildConflict(Prob,mode)
%
%  Prob   TOMLAB problem structure, describing a LP/QP/MILP/MIQP problem
%
%  mode   String indicating which type of conflict group set is desired.
%
%         A 'full' conflict group set will consist of one group for each
%         individual variable (upper+lower bound), linear, quadratic, sos
%         and logical constraint in the problem. This will be a very
%         large group set.
%
%         A 'minimal' set consists of at the most 6 groups: one each for all
%         variable lower+upper bounds, linear, sos, logical, quadratic
%         constraints.
%
% Inputs for (2): function confgrps = cpxBuildConflict(n,m_lin,m_quad,m_sos,m_log,'mode')
%
% n       Number of variables
% m_lin   Number of linear constraints
% m_quad  Number of quadratic constraints
% m_sos   Number of SOS constraints
% m_log   Number of logical constraints
% mode    Mode indicator as described above.
%

% Anders Goran, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2007 by Tomlab Optimization Inc., $Release: 11.0.0$
% Written Jan 20, 2006  Last modified Feb 3, 2006

function confgrps = cpxBuildConflict(P1,P2,P3,P4,P5,P6)

if nargin<6
   P6=[];
   if nargin<5
      P5=[];
      if nargin<4
         P4=[];
         if nargin<3
            P3=[];
            if nargin<2
               P2=[];
               if nargin < 1
                  error('TOMLAB:CPLEX:INPUT','cpxBuildConflict need the Prob structure as input');
               end
            end
         end
      end
   end
end

if(nargin<=2)
   %% Assume call syntax (1)
   if ~isstruct(P1)
      error('TOMLAB:CPLEX:INPUT','cpxBuildConflict need the Prob structure as input');
   end

   Prob = P1;
   mode = P2;

   nvars = Prob.N;
   m_lin = Prob.mLin;
   m_qc  = 0;
   m_log = 0;
   m_sos = 0;

   % Quadratic constraints in Prob.QP.qc
   qc = DefPar(Prob.QP,'qc');
   if ~isempty(qc), m_qc = length(qc); end

   % Fields in Prob.CPLEX
   CPLEX = DefPar(Prob,'CPLEX');
   ind = DefPar(CPLEX,'logcon');
   if ~isempty(ind), m_log = length(ind); end

   % Fields in Prob.MIP
   MIP = DefPar(Prob,'MIP');
   sos1 = DefPar(MIP,'sos1');
   sos2 = DefPar(MIP,'sos2');

   m_sos = length(sos1) + length(sos2);
else
   % Call syntax (2)
   nvars  = P1;
   m_lin  = P2;
   m_qc   = P3;
   m_sos  = P4;
   m_log  = P5;
   mode   = P6;
end

if isempty(mode), mode = 'full'; end

% Total number of conflict member entities in a "full" conflict refinement
n_tot = 2*nvars + m_lin + m_qc + m_sos + m_log;

switch(mode)
   case 'full',
      ngrps = n_tot;

      i = 1;
      for k=1:nvars
         confgrps(i).lowercol = k;
         i=i+1;
      end

      for k=1:nvars
         confgrps(i).uppercol = k;
         i=i+1;
      end

      for k=1:m_lin
         confgrps(i).linear = k;
         i=i+1;
      end

      for k=1:m_qc
         confgrps(i).quad = k;
         i=i+1;
      end

      for k=1:m_log
         confgrps(i).logical = k;
         i=i+1;
      end

      for k=1:m_sos
         confgrps(i).sos = k;
         i=i+1;
      end

   case 'minimal',
      confgrps(1).lowercol = 1:nvars;
      confgrps(2).uppercol = 1:nvars;

      if m_lin > 0
         confgrps(end+1).linear = 1:m_lin;
      end

      if m_qc > 0
         confgrps(end+1).quad = 1:m_qc;
      end

      if m_sos > 0
         confgrps(end+1).sos = 1:m_sos;
      end

      if m_log
         confgrps(end+1).indicator = 1:m_log;
      end
end






