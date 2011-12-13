%% LogMIP User's Manual Example 1a - Job Scheduling
% TomSym implementation of GAMS Example (LOGMIP1A,SEQ=332)
%
% Three jobs (A,B,C) must be executed sequentially in three steps, but not
% all jobs require all the stages. The objective is to obtain the sequence
% of tasks which minimizes the completion time. Once a job has started it
% cannot be interrupted. The objective is to obtain the sequence of task,
% which minimizes the completion time.
%
% Ref: Raman & Grossmann, Comp. & Chem. Eng., 18, 7, p.563-578, 1994.

toms 6x1 I
toms 3x1 J X

% Binary variables
toms 6x1 integer Y
cbnd1 = {0 <= Y <= 1};

toms T

% X and T are positive
cbnd2 = {0 <= X; 0 <= T};

eq1 = {T >= X(1) + 8
    T >= X(2) + 5
    T >= X(3) + 6};

eq2 = {X(1)-X(3) <= -5
    X(3)-X(1) <= -2
    X(2)-X(3) <= -1
    X(3)-X(2) <= -6
    X(1)-X(2) <= -5
    X(2)-X(1) <= 0};

objective = T;

cbnd3 = {X <= 20};

% Disjunction
eq3 = {(X(1)-X(3))*Y(1) <= -5*Y(1)
    (X(3)-X(1))*Y(2) <= -2*Y(2)
    (X(2)-X(3))*Y(3) <= -1*Y(3)
    (X(3)-X(2))*Y(4) <= -6*Y(4)
    (X(1)-X(2))*Y(5) <= -5*Y(5)
    (X(2)-X(1))*Y(6) <= 0*Y(6)};

eq4 = {Y(1) + Y(2) == 1
    Y(3) + Y(4) == 1
    Y(5) + Y(6) == 1};

options = struct;
options.solver = 'minlpBB';
constr = {cbnd1;cbnd2;cbnd3; eq1;eq3;eq4};
solution = ezsolve(objective,constr,[],options);