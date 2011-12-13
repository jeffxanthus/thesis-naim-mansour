% Limited demo testing
%
% inputs: n    # of variables
%         ncon # of constraints

function cpxDemo(n,ncon)

c=ones(n,1);

% Create ncon constraints, most of them redundant
A=ones(ncon,n);
b_L=zeros(ncon,1);
b_U=n*(1:ncon);

x_L=zeros(n,1);
x_U=ones(n,1); 

[x,slack,v,rc,f,ninf,sinf,Inform]=cplex(c,A,x_L,x_U,b_L,b_U); 

f,Inform