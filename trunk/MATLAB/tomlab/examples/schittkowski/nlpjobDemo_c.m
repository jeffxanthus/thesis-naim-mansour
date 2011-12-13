function cx = nlpjobDemo_c(x,Prob)

if Prob.USER.P == 1
   cx(1) = 9 - x(1)^2 - x(2)^2;
   cx(2) = 1 - x(1) - x(2);
else
end