function fx = nlpjobDemo_f(x,Prob)

if Prob.USER.P == 1
   fx(1) = (x(1) + 3)^2 + 1;
   fx(2) = x(2);
else
end