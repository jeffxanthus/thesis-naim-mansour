function J = nlpinfDemo_J(x,Prob)

if Prob.USER.P == 1
   J(1,1) = 2*x(1) + x(2);
   J(1,2) = 2*x(2) + x(1);
   J(2,1) = cos(x(1));
   J(2,2) = 0.0;
   J(3,1) = 0.0;
   J(3,2) = sin(x(2));
end
