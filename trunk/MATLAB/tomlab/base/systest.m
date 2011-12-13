% systest.m
%
% function systest(solvTypes,PriLevOpt,PriLev,wait);
%
% systest runs big test to check for bugs in TOMLAB
%
% INPUT:
% solvTypes   Which solvTypes to test
%             1=uc, 2=qp, 3=con, 4=ls, 5=lls, 6=cls, 7=mip, 8=lp,
%             9=glb, 10=glc
%
% solvType = 3 will first test conSolve, then nlpSolve
% solvType = 6 will first test cls_prob, then exp_prob
% If v3.0 then
%    test on minos for solvType 1 2 3 6 8
%    test on qpopt for solvType 2 8
%
% PriLevOpt   Printing level in the Solver. Default 2, short info each iter.
% wait        wait=1 if pause after each problem. Default = 1.
% PriLev      Printing level in PrintResult. Default full information == 5.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Nov 2, 1998.    Last modified June 4, 2005.

function systest(solvTypes,PriLevOpt,PriLev,wait)

if nargin < 4
   wait=[];
   if nargin < 3
      PriLev=[];
      if nargin < 2
         PriLevOpt=[];
         if nargin < 1
            solvTypes = [1 2 3 4 5 6 7 8 9 10];
         end
      end
   end
end

if isempty(PriLev),    PriLev=1;            end
if isempty(wait),      wait=0;              end
if isempty(PriLevOpt), PriLevOpt=0;         end

[TomV, os, TV] = tomlabVersion;

for i=1:length(solvTypes)
    solvType=solvTypes(i);
    if solvType == 1
       if TV(1)
          runtest('ucSolve',0,'uc_prob',[],PriLevOpt,wait,PriLev);
          runtest('sTrustr',0,'uc_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(2)
          runtest('minos',0,'uc_prob',[],PriLevOpt,wait,PriLev);
       end
    elseif solvType == 2
       if TV(1)
          runtest('qpSolve',0,'qp_prob',[],PriLevOpt,wait,PriLev);
          runtest('qld',0,'qp_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(2)
          runtest('qp-minos',0,'qp_prob',[],PriLevOpt,wait,PriLev);
          runtest('qpopt',0,'qp_prob',[],PriLevOpt,wait,PriLev);
          runtest('minos',0,'qp_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(4)
          runtest('snopt',0,'qp_prob',[],PriLevOpt,wait,PriLev);
          runtest('sqopt',0,'qp_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(3)
          runtest('npsol',0,'qp_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(6)
          runtest('bqpd',0,'qp_prob',[],PriLevOpt,wait,PriLev);
       end
    elseif solvType == 3
       if TV(1)
          runtest('conSolve',0,'con_prob',[],PriLevOpt,wait,PriLev);
          runtest('nlpSolve',0,'con_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(2)
          runtest('minos',0,'con_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(4)
          runtest('snopt',0,'con_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(3)
          runtest('npsol',0,'con_prob',[],PriLevOpt,wait,PriLev);
       end
    elseif solvType == 4
       if TV(1)
          runtest('clsSolve',1,'ls_prob',[],PriLevOpt,wait,PriLev);
          runtest('clsSolve',4,'ls_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(3)
          runtest('nlssol',0,'ls_prob',[],PriLevOpt,wait,PriLev);
       end
    elseif solvType == 5
       if TV(1)
          runtest('lsei',0,'lls_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(3)
          runtest('lssol',0,'lls_prob',[],PriLevOpt,wait,PriLev);
       end
    elseif solvType == 6
       if TV(1)
          runtest('clsSolve',1,'cls_prob',[],PriLevOpt,wait,PriLev);
          runtest('clsSolve',1,'exp_prob',[],PriLevOpt,wait,PriLev);
          runtest('conSolve',0,'cls_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(2)
          runtest('minos',0,'cls_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(4)
          runtest('snopt',0,'cls_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(3)
          runtest('nlssol',0,'exp_prob',[],PriLevOpt,wait,PriLev);
       end
    elseif solvType == 7
       Prob.optParam.IterPrint=0;
       if TV(1)
          runtest('mipSolve',0,'mip_prob',1:16,PriLevOpt,wait,PriLev,Prob);
       end
    elseif solvType == 8
       Prob.optParam.IterPrint=0;
       if TV(1)
          runtest('lpSimplex',0,'lp_prob',[],PriLevOpt,wait,PriLev,Prob);
          runtest('qld',0,'lp_prob',[],PriLevOpt,wait,PriLev,Prob);
       end
       if TV(2)
          runtest('minos',0,'lp_prob',[],PriLevOpt,wait,PriLev);
          runtest('lp-minos',0,'lp_prob',[],PriLevOpt,wait,PriLev);
          runtest('qpopt',0,'lp_prob',[],PriLevOpt,wait,PriLev);
       end
       if TV(4)
          runtest('snopt',0,'lp_prob',[],PriLevOpt,wait,PriLev);
          runtest('sqopt',0,'lp_prob',[],PriLevOpt,wait,PriLev);
       end
    elseif solvType == 9
       %runtest('glbSolve',0,'glb_prob',[],PriLevOpt,wait,PriLev);
       if TV(1)
          runtest('glbFast',0,'glb_prob',[1:32],PriLevOpt,wait,PriLev);
          runtest('glcFast',0,'glb_prob',[1:32],PriLevOpt,wait,PriLev);
          runtest('glbSolve',0,'glb_prob',[1:9,13:28],PriLevOpt,wait,PriLev);
          runtest('glcSolve',0,'glb_prob',[1:9,13:28],PriLevOpt,wait,PriLev);
          runtest('glcCluster',0,'glb_prob',[1:32],PriLevOpt,wait,PriLev);
       end
    elseif solvType == 10
       if TV(1)
          runtest('glcFast',0,'glc_prob',[1:31],PriLevOpt,wait,PriLev);
          runtest('glcSolve',0,'glc_prob',[1:31],PriLevOpt,wait,PriLev);
          runtest('glcCluster',0,'glc_prob',[1:22,24:31],PriLevOpt,wait,PriLev);
       end
    end
end
     
% MODIFICATION LOG
%
% 981102 hkh  Written
% 990616 hkh  Added constrained global optimization
% 990704 hkh  Added linear programming
% 020105 hkh  Adding new global optimization solvers
% 020110 hkh  Must not run conSolve on exp fitting with analytic Hessian!
% 020701 hkh  Use new way of checking Tomlab solver options
% 020820 ago  Error if TV(5)==1 - this is CGO, not BQPD
% 040402 hkh  Avoid problem 23 in glc_prob for glcCluster, pure IP  
% 050117 med  mlint revision
% 050604 hkh  Change lpSolve to lpSimplex
