% function Prob = establishingacollegetimetable(lessons1, lessons2,
% subject, slots)
%
% Creates a TOMLAB MILP problem for establishing a college time table
%
% INPUT PARAMETERS
%
% OUTPUT PARAMETERS
% Prob          A TOMLAB problem defining the problem (type MIP).

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Dec 2, 2005.   Last modified Dec 2, 2005.

function Prob = establishingacollegetimetable(lessons1, lessons2, subject, slots)

if nargin < 4
   error('The function requires 4 inputs');
end

if isempty(lessons1) | isempty(lessons2) | isempty(subject) | isempty(slots)
   error('One of the inputs are empty');
end

n1    = length(lessons1);        %teachers
n2    = 2;                       %2
n3    = slots;                   %slots

n     = n1*n2*n3; %teachers (class 1), teachers (class 2).... for all slots

% FORMULATE PROBLEM
% All variables are binary
IntVars   = ones(n,1);
x_L       = zeros(n,1);
x_U       = ones(n,1);

% Course constraints
b_L1 = [lessons1;lessons2];
b_U1 = b_L1;
A1   = zeros(n1*n2,n);
for i=1:n1
   A1(i,i:n1*n2:n-n1*n2+i) = ones(1,n3);
end
for i=1:n1
   A1(n1+i,i+n1:n1*n2:n-n1+i) = ones(1,n3);
end

% Teacher constraint, one teacher per slot and class
b_L2 = -inf*ones(n2*n3,1);
b_U2 = ones(n2*n3,1);
A2   = zeros(n2*n3,n);
for i=1:n2*n3
   A2(i,(i-1)*n1+1:i*n1) = ones(1,n1);
end

% Class constraint, one class per slot
b_L3 = -inf*ones(n1*n3,1);
b_U3 = ones(n1*n3,1);
A3   = zeros(n1*n3,n);
for i=1:n1
   for j=1:n3
      A3(j+(i-1)*n3,[(j-1)*n1*n2+i, (j-1)*n1*n2+i+n1]) = ones(1,2);
   end
end

% At least one 2-hour slot of a subject per day

% Biology
A4 = zeros(5*2,n);
counter = 1;
for i=1:5
   ind = 2;
   A4(counter, [(i-1)*n1*n2*4+ind:n1*n2:i*n1*n2*4]) = ones(1,4);
   counter = counter + 1;
   A4(counter, [(i-1)*n1*n2*4+ind+n1:n1*n2:i*n1*n2*4+n1]) = ones(1,4);
   counter = counter + 1;
end
b_L4 = -inf*ones(10,1);
b_U4 = ones(10,1);

% History
A5 = zeros(5*2,n);
counter = 1;
for i=1:5
   ind = 3;
   A5(counter, [(i-1)*n1*n2*4+ind:n1*n2:i*n1*n2*4]) = ones(1,4);
   counter = counter + 1;
   A5(counter, [(i-1)*n1*n2*4+ind+n1:n1*n2:i*n1*n2*4+n1]) = ones(1,4);
   counter = counter + 1;
end
b_L5 = -inf*ones(10,1);
b_U5 = ones(10,1);

% Physics
A6 = zeros(5*2,n);
counter = 1;
for i=1:5
   ind = 6;
   A6(counter, [(i-1)*n1*n2*4+ind:n1*n2:i*n1*n2*4]) = ones(1,4);
   counter = counter + 1;
   A6(counter, [(i-1)*n1*n2*4+ind+n1:n1*n2:i*n1*n2*4+n1]) = ones(1,4);
   counter = counter + 1;
end
b_L6 = -inf*ones(10,1);
b_U6 = ones(10,1);

% Math 1
A7 = zeros(5,n);
counter = 1;
for i=1:5
   ind = 5;
   A7(counter, [(i-1)*n1*n2*4+ind:n1*n2:i*n1*n2*4]) = ones(1,4);
   counter = counter + 1;
end
b_L7 = -inf*ones(5,1);
b_U7 = ones(5,1);

% Math 2
A8 = zeros(5,n);
counter = 1;
for i=1:5
   ind = 4;
   A8(counter, [(i-1)*n1*n2*4+ind+n1:n1*n2:i*n1*n2*4+n1]) = ones(1,4);
   counter = counter + 1;
end
b_L8 = -inf*ones(5,1);
b_U8 = ones(5,1);

% Sport has to be on thursday afternoon
x_L(8+n1*n2*14,1) = 1;
x_L(9+n1+n1*n2*14,1) = 1;

% No lessons Monday morning
x_U(1:n1*n2,1) = zeros(n1*n2,1);

% Teacher 4 is away on Monday morning
x_U([n1+4,n1+4+n1*n2],1) = zeros(2,1);

% Teacher 2 is away on Wednesdays
x_U(8*n1*n2+2:n1:11*n1*n2+2+n1,1) = zeros(8,1);

% Merge constraints
A   = [A1;A2;A3;A4;A5;A6;A7;A8];
b_L = [b_L1;b_L2;b_L3;b_L4;b_L5;b_L6;b_L7;b_L8];
b_U = [b_U1;b_U2;b_U3;b_U4;b_U5;b_U6;b_U7;b_U8];

reper = [ones(n1*n2,1); zeros(n1*n2*2,1); ones(n1*n2,1)];
c   = repmat(reper, 5 , 1);
Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], 'Establishing a College Time Table', [], [], IntVars);

% MODIFICATION LOG
%
% 051202 med   Created
% 060123 med   Added constraints for days (one subject per day)