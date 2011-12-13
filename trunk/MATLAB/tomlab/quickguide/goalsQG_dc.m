% function dc=goalsQG_dc(x, Prob)
%
% The goals problems are Multi Criterium Unconstrained & Constrained Nonlinear Problems
% 'EASY-FIT TP355'
%
% goalsQG_dc computes the gradient to the nonlinear constraints c in the
% point x.

function dc=goalsQG_dc(x, Prob)

dc1=30*x(2)^2*x(4)^2*x(3)+60*x(2)*x(4)*x(1)-12*x(2)*x(4)*x(3)-126*x(2)^2*x(4)^2*x(1)-6*x(1)...
    +60*x(4)-60*x(2)-30*x(1)*x(4)^2+6*x(3)*x(4)^2+300*x(2)^2*x(4);
dc2=600*x(2)*x(4)*x(1)-6*x(2)*x(4)^2*x(3)^2-12*x(1)*x(3)*x(4)-120*x(2)*x(4)*x(3)-126*x(2)*x(4)^2*x(1)^2 ...
    -60*x(1)+6*x(4)-600*x(2)-30*x(2)*x(4)^2+30*x(1)^2*x(4)+60*x(2)*x(4)^2*x(3)*x(1);
dc3=-6*x(2)^2*x(4)^2*x(3)-12*x(2)*x(4)*x(1)+30*x(2)^2*x(4)^2*x(1)+1.5*x(3)+6*x(1)*x(4)^2-60*x(2)^2*x(4);
dc4=-6*x(2)^2*x(4)*x(3)^2+12*x(1)*x(3)*x(4)-126*x(2)^2*x(4)*x(1)^2-12*x(1)*x(2)*x(3)+60*x(1)-6*x(4)+6*x(2)...
    -30*x(2)^2*x(4)-30*x(1)^2*x(4)+60*x(2)^2*x(4)*x(3)*x(1)+30*x(1)^2*x(2)-60*x(2)^2*x(3)+300*x(1)*x(2)^2;
dc=[dc1 dc2 dc3 dc4];