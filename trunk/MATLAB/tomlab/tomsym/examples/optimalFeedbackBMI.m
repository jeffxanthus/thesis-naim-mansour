%% Optimal Feedback - A tomSym BMI demonstration
% A sample problem from the PENBMI manual.
%
% LQ optimal feedback
%
% Minimize
%   trace(P)
% subject to
%   P >= 0                                  (LMI)
%   (A+B*K)'*P+P*(A+B*K) <= -eye(2)-K'*K    (BMI)

% Variables
toms 2x2 symmetric P
toms 1x2           K

% Constants
A = [-1 2;-3 -4];
B = [1;1];

% Constraints
F = { MI( P >= 0 )
    MI( (A+B*K)'*P + P'*(A+B*K) <= -eye(2) - K'*K )};

%Solve the problem
solution = ezsolve(trace(P),F,[],'Optimal feedback');

% Evaluate K using the returned solution
K_opt = subs(K,solution)