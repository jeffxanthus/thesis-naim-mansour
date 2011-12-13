% nllsQG_J - Jacobian matrix for nnls quick guide
%
% function J = nllsQG_J(x, Prob)

function J = nllsQG_J(x, Prob)

y=Prob.LS.y;
t=Prob.LS.t;
m=size(y,1);
uP = Prob.uP;
a=uP(1)*x(1)/(x(3)*(x(1)-x(2)));
b=x(1)-x(2);
J=zeros(m,3*size(x,2));
for i=1:m
    e1=exp(-x(1)*t(i)); e2=exp(-x(2)*t(i));
    J(i,1)=a*(t(i)*e1+(e2-e1)*(1-1/b));
    J(i,2)=a*(-t(i)*e2+(e2-e1)/b);
    J(i,3)=-a*(e2-e1)/x(3);
end