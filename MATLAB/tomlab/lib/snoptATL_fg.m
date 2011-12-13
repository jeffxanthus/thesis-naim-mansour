% snoptA callback for traditional TOMLAB problems.
% To demonstrate (however inefficiently) how snoptA could
% interface with TOMLAB.

function [F,G,mode] = snoptATL_fg(x,status,Prob)

mode  = status(1);
needf = status(2);
needG = status(3);

if needf
   % Objective function
   F = nlp_f(x,Prob);
   % Nonlinear constraints
   if Prob.mNonLin>0
      F = [F ; nlp_c(x,Prob) ];
   end
   % Entirely linear functions have zero nonlinear contributions
   F(end+1:end+Prob.mLin)=0;
else
   F=[];
end

if needG
   % Objective gradient
   G = nlp_g(x,Prob)';
   % Nonlinear constraints Jacobian
   if Prob.mNonLin>0
      G = [G ; nlp_dc(x,Prob) ];
   end
   % Entirely linear functions have zero nonlinear jacobian
   G(end+1:end+Prob.mLin,end)=0;
else
   G=[];
end
