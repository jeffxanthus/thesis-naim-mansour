% function Result0 = runtest(Solver, SolverAlg, probFile, probNumbs,...
%                           PriLevOpt, wait, PriLev, Prob0)
% 
% runtest runs all selected problems defined in a problem file
% for the given solver
%
% INPUT: 
% Solver      Name of solver, default conSolve.
% SolverAlg   A vector of numbers defining which of the Solver algorithms to
%             try. For each element in SolverAlg, all probNumbs are solved
%             Leave empty, or set 0 if to use the default algorithm.
% probFile    Problem definition file. probFile is by default set to
%             con_prob if Solver is conSolve, uc_prob if Solver is ucSolve
%             and so on. 
% probNumbs   A vector with problem numbers to run. If empty, run all in file.
% PriLevOpt   Printing level in the Solver. Default 2, short info each iter.
% wait        wait=1 if pause after each problem. Default = 1.
% PriLev      Printing level in PrintResult. Default full information == 5.
% Prob0       Initial problem structure for each problem
%
% OUTPUT: 
% Result0     Structure matrix with some of the Result output

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Oct 28, 1998.  Last modified Nov 23, 2004.

function Result0 = runtest(Solver, SolverAlg, probFile, probNumbs,...
                          PriLevOpt, wait, PriLev, Prob0)

if nargin < 8
   Prob0 = [];
   if nargin < 7
      PriLev=[];
      if nargin < 6
         wait=[];
         if nargin < 5
            PriLevOpt=[];
            if nargin < 4
               probNumbs=[];
               if nargin < 3
                  probFile=[];
                  if nargin < 2
                     SolverAlg=[];
                     if nargin < 1
                        Solver=[];
end, end, end, end, end, end, end, end

if isempty(PriLev),    PriLev=2;            end
if isempty(wait),      wait=1;              end
if isempty(PriLevOpt), PriLevOpt=2;         end
if isempty(SolverAlg), SolverAlg=0;         end
if isempty(Solver),    Solver='conSolve';   end

% Find default file
if isempty(probFile)  
   if strcmpi(Solver,'conSolve')
      probFile = 'con_prob';
   else
      [SolvList,SolvTypeList] =  SolverList([],0,1);

      i=strmatch(lower(Solver),lower(SolvList),'exact');

      ss=str2mat('uc_prob','qp_prob','con_prob','ls_prob', 'lls_prob',...
            'cls_prob','mip_prob','lp_prob','glb_prob','glc_prob',...
            'miqp_prob','minlp_prob','sdp_prob','bmi_prob','exp_prob');
      probFile=deblank(ss(SolvTypeList(i),:));
   end
end

probNames = feval(probFile); % Names of available problems

if isempty(probNumbs)
   probNumbs=[1:size(probNames,1)];
end

ask=-1;

for i = 1 : length(SolverAlg);
    alg = SolverAlg(i);

    if PriLev == 1 & PriLevOpt <=0
       fprintf('Solver: %s',Solver);
       fprintf('. Algorithm %d\n',alg);
    end
   
    for j = 1 : length(probNumbs)
        P=probNumbs(j);
        if ~isempty(Prob0)
           Prob0.P=P;
           Prob0.probFile=[];
        end
        Prob = probInit(probFile,P,ask,Prob0);
        Prob.Solver.Alg=alg;

        if PriLev == 0 
           fprintf('Solver: %s',Solver);
           fprintf('. Algorithm %d\n',alg);
           fprintf('Problem %d - %s',P,probNames(P,:));
           fprintf('\n');
        end

        Prob.PriLevOpt=PriLevOpt;
        %if isfield(Prob.LS,'SepAlg')
        %   if Prob.LS.SepAlg==1
        %      n=round(length(Prob.x_0)/2);
        %      Prob.x_0=Prob.x_0(1:n);
        %      if ~isempty(Prob.x_L), Prob.x_L=Prob.x_L(1:n); end
        %      if ~isempty(Prob.x_U), Prob.x_U=Prob.x_U(1:n); end
        %      if ~isempty(Prob.A),   Prob.A=Prob.A(:,1:n);   end
        %   end
        %end

        Result = tomRun(Solver, probFile, P, Prob, PriLev, ask);
      
        if ~isempty(Result)
           if isnan(any(Result.x_k)) | isnan(Result.f_k)
              fprintf('\n *** WARNING NaN in output ***\n');
           end
           if isinf(any(Result.x_k)) | isinf(Result.f_k)
              fprintf('\n *** WARNING Inf in output ***\n');
           end
           if wait & length(probNumbs) > 1
              pause
           end
           Result0(i,j).Prob     =Result.Prob;
           Result0(i,j).x_k      =Result.x_k;
           Result0(i,j).f_k      =Result.f_k;
           Result0(i,j).FuncEv   =Result.FuncEv;
           Result0(i,j).GradEv   =Result.GradEv;
           Result0(i,j).ConstrEv =Result.ConstrEv;
           Result0(i,j).ResEv    =Result.ResEv;
           Result0(i,j).JacEv    =Result.JacEv;
           Result0(i,j).CPUtime  =Result.CPUtime;
           Result0(i,j).REALtime =Result.REALtime;
           Result0(i,j).Iter     =Result.Iter;
           Result0(i,j).p_dx     =Result.p_dx;
        else
           Result0(i,j).x_k      =[];
        end
    end
end

% MODIFICATION LOG
%
% 981109   hkh   Change default PriLev=2
% 981111   hkh   Update with LP solvers
% 981127   hkh   Added global optimization solver ego to list
% 981130   hkh   Use SolverList instead of hard coded values.
%                Generalized to use the correct default ???_prob file.
% 981202   mbk   Changes in comments.
% 990928   hkh   Added mipRun
% 991010   hkh   Added initial structure Prob0. Handle separable NLLS
% 000830   hkh   New driver format in v3.0
% 001009   hkh   Removed saving of flops
% 001023   hkh   Use tomRun as driver
% 040102   hkh   Mideva part in comments, should not be needed any longer
% 040126   hkh   Changed handling of empty arguments, or probFile empty
% 040126   hkh   Added probTypes 11-14
% 041123   hkh   Change call to tomRun

