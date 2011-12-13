% snoptA_Ex1_fg
%
% User function callback for snoptA_Ex1 example
%
% Calculates the nonlinear components f(x) of the function 
% vector F(x)=f(x)+A*x and also the Jacobian of f(x)  
%

function [f,G,mode] = snoptA_Ex1_fg(x,status,Prob)
  
mode  = status(1); % 1=first call
needf = status(2);
needG = status(3);

% Function vector
if needf>0
   f = [ ...
      exp(x(2))*x(4)+x(2)^2 ; ...
      x(3)^2+sin(x(4)) ; ...
      0 ...
      ];
else
   f = [];
end

% Nonlinear Jacobian (including a row for the objective gradient)
if needG>0
   i = [1,2,1,2];
   j = [2,3,4,4];
   v = [exp(x(2))*x(4)+2*x(2),2*x(3),exp(x(2)),cos(x(4))];
   G = sparse(i,j,v,3,5);
   % Does not have to be sparse, can even change sparse/dense during a
   % single run!
else
   G = [];
end
