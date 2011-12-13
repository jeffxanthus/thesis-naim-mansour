probFile = 'ucnew_prob';           % Problem definition file.
P        = 18;                     % Problem number.
Prob     = probInit(probFile, P);  % Setup Prob structure.
Result   = ucSolve(Prob);
PrintResult(Result);