function J = nlpinfDemo_J(x,Prob)

if Prob.USER.P == 1
   J(1,1) = -20*x(1);
   J(1,2) = 10;
   J(1,3) = 0;
   J(1,4) = 0;
   J(2,1) = -1;
   J(2,2) = 0;
   J(2,3) = 0;
   J(2,4) = 0;
elseif Prob.USER.P == 2
   t = Prob.LS.t;
   J(:,1) = t.*(t + x(2))./(t.^2 + x(3)*t + x(4));
   J(:,2) = x(1)*t./(t.^2 + x(3)*t + x(4));
   J(:,3) = -x(1)*(t.^2).*(t + x(2))./(t.^2 + x(3)*t + x(4)).^2;
   J(:,4) = -x(1)*t.*(t + x(2))./(t.^2 + x(3)*t + x(4)).^2;
end
