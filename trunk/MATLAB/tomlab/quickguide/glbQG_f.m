% glbQG_f - function value for Shekel 5, Quickguide
%
% function f = glbQG_f(x, Prob)

function f = glbQG_f(x, Prob)

A = [4 4 4 4;1 1 1 1;8 8 8 8;6 6 6 6;3 7 3 7];
c = [.1 .2 .2 .4 .4]';
f=0;
for i = 1:5
    f = f - 1./( (x-A(i,:)' )'*( x-A(i,:)' ) + c(i) );
end