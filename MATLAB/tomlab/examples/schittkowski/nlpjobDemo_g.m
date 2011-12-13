function gx = nlpjobDemo_g(x,Prob)

if Prob.USER.P == 1
   gx(1,1) = 2*(x(1) + 3);
   gx(1,2) = 0;
   gx(2,1) = 0;
   gx(2,2) = 1;
else
end