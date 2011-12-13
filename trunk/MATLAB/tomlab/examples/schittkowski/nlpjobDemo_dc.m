function dcx = nlpjobDemo_dc(x,Prob)

if Prob.USER.P == 1
   dcx(1,1) = -2*x(1);
   dcx(1,2) = -2*x(2);
   dcx(2,1) = -1;
   dcx(2,2) = -1;
else
end
