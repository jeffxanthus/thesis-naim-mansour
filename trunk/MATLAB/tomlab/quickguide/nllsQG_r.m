% nllsQG_r - function value for nlls quick guide
%
% function r = nllsQG_r(x, Prob)

function r = nllsQG_r(x, Prob)

y = Prob.LS.y;
t = Prob.LS.t;
m = size(y,1);
uP = Prob.uP;
r=zeros(m,size(x,2));
for i=1:m
   r(i)=uP(1)*x(1)*(exp(-x(2)*t(i))-exp(-x(1)*t(i)))/(x(3)*(x(1)-x(2)));
end
r=r-y;