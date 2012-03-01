function [ x ] = BasisPursuitFast( A, b, Bcon, Ccon, thetaClipPos, thetaClipNeg )
%BasisPursuitFast Regular basis pursuit algorithm. Solves LP
%   Detailed explanation goes here


%%%%%%%%%%%%%%%%%
% Dees werkt, maar kan mss nog iets beter...
%%%%%%%%%%%%%%%%
[m,n] = size(A); 
Amatrix = [-Ccon' Bcon' A' 2*eye(n)]; 
bvector = ones(n,1); 
cvector = [-thetaClipNeg; -thetaClipPos; -b; zeros(n,1)]; 
lb = [zeros(m,1); zeros(m,1); -inf*ones(m,1); zeros(n,1)]; 
ub = [inf*ones(m,1); inf*ones(m,1); inf*ones(m,1); ones(n,1)];


%Acon = [];
%bcon = [];

% thetaClipPos
% thetaClipNeg
% 
% size(Amatrix)
% size(bvector)
% size(Acon)
% size(bcon)
% size(cvector)

% track the time required to solve the dual problem profile on 
[primal,obj,exitflag,output,dual] = linprog(cvector,[],[],Amatrix,bvector,lb,ub); 

% x is the dual variable corresponding to the equality constraints in dual problem 
x = dual.eqlin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% [m,n] = size(A); 
% Amatrix = [-Ccon' Bcon' A' 2*eye(n)]; 
% bvector = ones(n,1); 
% cvector = [-thetaClipNeg; thetaClipPos; -b; zeros(n,1)]; 
% lb = [zeros(m,1); zeros(m,1); -inf*ones(m,1); zeros(n,1)]; 
% ub = [inf*ones(m,1); inf*ones(m,1); inf*ones(m,1); ones(n,1)];
% Acon = [Ccon' -Bcon' zeros(n) zeros(n)];
% bcon = bvector;
% 
% %Acon = [];
% %bcon = [];
% 
% % thetaClipPos
% % thetaClipNeg
% % 
% % size(Amatrix)
% % size(bvector)
% % size(Acon)
% % size(bcon)
% % size(cvector)
% 
% % track the time required to solve the dual problem profile on 
% [primal,obj,exitflag,output,dual] = linprog(cvector,Amatrix,bvector,[],[],lb,ub); 
% 
% % x is the dual variable corresponding to the equality constraints in dual problem 
% x = dual.ineqlin;
% size(x)
% n
% x = x(1:n);


end

