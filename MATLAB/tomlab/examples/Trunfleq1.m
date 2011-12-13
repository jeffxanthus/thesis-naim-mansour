% function [fval, exitflag, output, x] = Trunfleq1
%
% The runfleq1 example in Optimization toolbox is using the
% option HessMult.
%
% This is not efficient in Tomlab, because no solver in Tomlab is utilizing
% the Hessian information in this way.
%
% Instead solve the problem with explicit Hessian

function [fval, exitflag, output, x] = Trunfleq1

% Same data as fleq1.mat in Math Works Optimization Toolbox

load('fleq1xx.mat')     % Get V, Aeq, beq, 
n = 1000;           % problem dimension

% Initial starting point
xstart = -ones(n,1); xstart(2:2:n,1) = ones(length(2:2:n),1); 

options = optimset('GradObj','on','Hessian','on', 'Display','iter'); 
if exist('tomlablic')
   global otxProb
   otxProb = ProbDef;
   otxProb.HessPattern=speye(1000,1000);
   for i=1:n-1
       otxProb.HessPattern(i,i+1) = 1;
       otxProb.HessPattern(i+1,i) = 1;
   end
   % otxProb.Solver.Tomlab = 'conopt';
   otxProb.Solver.Tomlab = 'knitro';
end
tic
[x,fval,exitflag,output] = fmincon(@brownfghxx, ...
   xstart,[],[],Aeq,beq,[],[], [],options,V);
toc