[TomV,os, TV] = tomlabVersion;
diary testtom.txt
if 1
   % Test of opt tb 2.0 compatible interfaces
   testfminsearch
   testfmincon
   testfminunc
   testlinprog
   testlsqcurvefit
   testlsqlin
   testlsqnonlin
   testquadprog
end
if 1
   testlsqnonneg
end

if 1
   % Test 1st problem for each problem type and Tomlab solver
   tomRun('ucSolve','uc_prob',1,2);
   tomRun('qpSolve','qp_prob',1,2);
   tomRun('conSolve','con_prob',1,2);
   tomRun('clsSolve','ls_prob',1,2);
   tomRun('clsSolve','cls_prob',1,2);
   tomRun('mipSolve','mip_prob',1,2);
   tomRun('lpSolve','lp_prob',1,2);
   tomRun('glbSolve','glb_prob',1,2);
   tomRun('glcSolve','glc_prob',1,2);
end

systest(1:4)

systest(5);


% Test exponential problems, avoiding problem 46-51. 
% The model is not implemented yet for these problems
%runtest('clsSolve',0,'exp_prob',1:45,1,0,1);

systest(6);

systest(7);

systest(8);

% Test glb problems, avoiding problem 10-12. Takes too long time
systest(9);

%runtest('glbSolve',0,'glb_prob',1:9,1,0,1);
%runtest('glbSolve',0,'glb_prob',13:28,1,0,1);

% Test global constrained problems
systest(10);


PAUS = 0;
disp('The next two demo files, democon and demouc, have pause statements');
democon
demouc

diary off
