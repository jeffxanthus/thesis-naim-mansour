% lpconQG_dc - nonlinear constraint gradient matrix for lpcon quick guide
%
% function dc = lpconQG_dc(x, Prob)

function dc = lpconQG_dc(x, Prob)

% Circle-triangle
h1=x(1)^2+x(2)^2;
h2=x(1)*x(3)+x(2)*x(4);
if h1==0
    dc=[1E10,1E10,1E10,1E10];
else
    dc=[2*h2*x(3)/h1-2*x(1)*(h2/h1)^2, 2*h2*x(4)/h1-2*x(2)*(h2/h1)^2,...
        2*h2*x(1)/h1-2*x(3), 2*h2*x(2)/h1-2*x(4)];
end