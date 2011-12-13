
lpExample;

% linprog needs linear inequalities and equalities to be given separately
% If the problem has both linear inequalities (only upper bounded) 
% and equalities we can easily detect which ones doing the following calls

ix = b_L==b_U;
E  = find(ix);
I  = find(~ix);

[x, fVal, ExitFlag, Out, Lambda] = linprog(c, A(I,:),b_U(I),...
    A(E,:), b_U(E), x_L, x_U, x_0);

% If the problem has linear inequalites with different lower and upper bounds
% the problem can be transformed using the TOMLAB routine cpTransf.
% See the example file tomlab\examples\testlinprog.m for an example.

fprintf('\n');
fprintf('\n');
disp('Run TOMLAB linprog on LP Example');
fprintf('\n');
xprinte(A*x-b_U,          'Constraints Ax-b_U ');
xprinte(Lambda.lower,     'Lambda.lower:      ');
xprinte(Lambda.upper,     'Lambda.upper:      ');
xprinte(Lambda.eqlin,     'Lambda.eqlin:      ');
xprinte(Lambda.ineqlin,   'Lambda.ineqlin:    ');
xprinte(x,                'x:                 ');
format compact
disp('Output Structure')
disp(Out)
fprintf('Function value %30.20f. ExitFlag %d\n',fVal,ExitFlag);

