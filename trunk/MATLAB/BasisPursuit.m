function [ x ] = BasisPursuit( A, b, Bcon, Ccon, thetaClipPos, thetaClipNeg )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Solves ?||x||1 + 1 (||Ax ? b||2)2 ?2
% where ? = ?	2log(p) % Basis pursuit with denoising
[m,n] = size(A); 
gamma = 0.003;
% Our QP is %Min??n (ui+vi)+1pTp
% subject to A(u ? v) ? Ip = b, u ? 0andv ? 0
% We split p = (p1 ? p2) in the actual formulation

H = zeros(2*m+2*n,2*m+2*n); 
H(1:2*m,1:2*m) = [eye(m) -eye(m); -eye(m) eye(m)];

f = [zeros(2*m,1); gamma*ones(2*n,1)];
% Aeq = [-eye(m) eye(m) A -A; eye(m) -eye(m) -A A; zeros(2*m,n) zeros(2*m,n) -eye(2*m,2*n)]; 
% beq = [b; -b; zeros(n,1); zeros(n,1)];

size(A)
size(Bcon)
size(Ccon)
sum(sum(abs(A)))
sum(sum(abs(Bcon)))
sum(sum(abs(Ccon)))

sum(sum(Bcon-Ccon))


Aeq = [-eye(m) eye(m) A -A; eye(m) -eye(m) -A A; zeros(2*m,n) zeros(2*m,n) -eye(2*m,2*n); zeros(m) zeros(m) Bcon -Bcon; zeros(m) zeros(m) -Ccon Ccon]; 
beq = [b; -b; zeros(n,1); zeros(n,1); thetaClipPos; thetaClipNeg];

% Aeq = [-eye(m) eye(m) A -A; eye(m) -eye(m) -A A]; 
% beq = [b; -b];

size(Aeq)
size(beq)

A = [zeros(m) zeros(m) -eye(m,2*n)];
B = zeros(m,1);

%sol = sqsolqp(H,Aeq,beq,f,1e-2);

%opt = optimset('TolPCG', 1e-2, 'TolX',1e-2, 'TolFun',1e-2, 'Algorithm', 'interior-point-convex');
%sol = quadprog(H,f,A,B,Aeq,beq,-inf,inf, [], opt); 
% x is (u-v)

x_L = ones(size(Aeq,2),1)*(-2000);
x_U = ones(size(Aeq,2),1)*1000;
b_U = ones(size(Aeq,1),1)*2000;

[sol, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, ...
          glnodes, confstat, iconfstat, sa, cpxControl, presolve] = ...
          cplex(f, Aeq, x_L, x_U, beq, b_U, [], [], [], [], [],[], [], [], [], [], H);
      

     
x = sol(2*m+1:2*m+n) - sol(2*m+n+1:2*m+2*n);

end

