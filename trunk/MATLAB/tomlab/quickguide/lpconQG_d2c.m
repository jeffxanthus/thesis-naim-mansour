% lpconQG_d2c - The second part of the Hessian to the Lagrangian function for the
%               nonlinear constraints for lpcon quick guide, i.e.:
%
%           lam' * d2c(x)
%
% in
%
%   L(x,lam) =   f(x) - lam' * c(x)
% d2L(x,lam) = d2f(x) - lam' * d2c(x) = H(x) - lam' * d2c(x)
% 
% function d2c = lpconQG_d2c(x, lam, Prob)

function d2c=lpconQG_d2c(x, lam, Prob)

% Circle-triangle
if lam(1) ~=0
    d2c=Inf*ones(4,4);
    if x(1)==0 & x(2)==0
        return
    end
    h1=x(1)^2+x(2)^2;
    h2=x(1)*x(3)+x(2)*x(4);
    d2c(1,1)=2*x(3)^2/h1-4*x(1)*x(3)*h2/h1^2-...
        2*(h2/h1)^2-4*x(1)*x(3)*h2/h1^2+8*x(1)^2*h2^2/h1^3;
    d2c(1,2)=2*x(3)*x(4)/h1-4*x(2)*x(3)*h2/h1^2-...
        2*(h2/h1)^2-4*x(1)*x(4)*h2/h1^2-8*x(1)*x(2)*h2^2/h1^3;
    d2c(1,3)=2*x(1)*x(3)/h1+2*h2/h1-4*x(1)^2*h2/h1^2;
    d2c(1,4)=2*x(2)*x(3)/h1-4*x(1)*x(2)*h2/h1^2;
    d2c(2,2)=2*x(4)^2/h1-4*x(2)*x(4)*h2/h1^2-...
        2*(h2/h1)^2-4*x(2)*x(4)*h2/h1^2+8*x(2)^2*h2^2/h1^3;
    d2c(2,3)=2*x(1)*x(4)/h1-4*x(1)*x(2)*h2/h1^2;
    d2c(2,4)=2*x(2)*x(4)/h1+2*h2/h1-4*x(2)^2*h2/h1^2;
    d2c(3,3)=2*x(1)^2/h1-2;
    d2c(4,4)=2*x(2)^2/h1-2;
    d2c(3,4)=2*x(1)*x(2)/h1;
    d2c(2,1)=d2c(1,2);
    d2c(3,1)=d2c(1,3);
    d2c(4,1)=d2c(1,4);
    d2c(3,2)=d2c(2,3);
    d2c(4,2)=d2c(2,4);
    d2c(4,3)=d2c(3,4);

    d2c=lam(1)*d2c;
else
    d2c=zeros(4,4);
end