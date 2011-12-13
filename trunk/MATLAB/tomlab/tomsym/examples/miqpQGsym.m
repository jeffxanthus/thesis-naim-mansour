% miqpQGsym - TomSym version of miqpQG

toms integer x
toms         y

objective = -6*x + 2*x^2 + 2*y^2 - 2*x*y;
constraints = {x+y<=1.9, x>=0, y>=0};

solution = ezsolve(objective,constraints);