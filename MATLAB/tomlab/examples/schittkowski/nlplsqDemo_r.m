function rx = nlpinfDemo_r(x,Prob)

if Prob.USER.P == 1
   rx(1) = 10*(x(2) - x(1)^2);
   rx(2) = 1 - x(1);
elseif Prob.USER.P == 2
   t = Prob.LS.t;
   y = Prob.LS.y;
   rx = [x(1)*t.*(t + x(2))./(t.^2 + x(3)*t + x(4)) - y];
end
