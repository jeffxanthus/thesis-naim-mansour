function rx = nlpinfDemo_r(x,Prob)

if Prob.USER.P == 1
   rx(1) = x(1)^2 + x(2)^2 + x(1)*x(2) - 1.0;
   rx(2) = sin(x(1));
   rx(3) = -cos(x(2));
end
