% function J = ls_J(x, Prob)
%
% Computes the Jacobian to least squares problem Gisela
%
% J(i,j) is dr_i/d_x_j

function J = ls1_J(x, Prob)

global US_A US_B


% Parameter K is input in the structure Prob

a  = Prob.userParam.K * x(1)/(x(3)*(x(1)-x(2)));
b  = x(1)-x(2);
t  = Prob.LS.t;

% Pick up the globally saved exponential computations
e1 = US_A;
e2 = US_B;

% Compute the three columns in the Jacobian, one for each of variable
J = a * [ t.*e1+(e2-e1)*(1-1/b), -t.*e2+(e2-e1)/b, (e1-e2)/x(3)];

% The loop version is given below, it is used in ls_prob.m, problem 3
% m = size(Prob.LS.y,1);
% J=zeros(m,3);
% for i=1:m
%     e1=exp(-x(1)*t(i)); e2=exp(-x(2)*t(i));
%     J(i,1)=a*(t(i)*e1+(e2-e1)*(1-1/b));
%     J(i,2)=a*(-t(i)*e2+(e2-e1)/b);
%     J(i,3)=-a*(e2-e1)/x(3);
% end