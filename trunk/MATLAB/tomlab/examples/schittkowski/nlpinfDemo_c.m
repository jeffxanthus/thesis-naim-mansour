function cx = nlpinfDemo_c(x, Prob)
if Prob.USER.P == 1
   cx = [];
elseif Prob.USER.P == 2
   t = Prob.LS.t;
   y = Prob.LS.y;
   cx = [x(1)*t(1)*(t(1) + x(2))/(t(1)^2 + x(3)*t(1) + x(4)) - y(1);
         x(1)*t(end)*(t(end) + x(2))/(t(end)^2 + x(3)*t(end) + x(4)) - y(end)];
end
