% Defines Hessian for mixed-integer nonlinear programming (MINLP) problems.
%
% function H = minlp_H(x,Prob)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function H = minlp_H(x,Prob)

x=x(:);
P = Prob.P;

if P==1 | P==10
    H = spalloc(5,5,0);
elseif P==2
    H = spalloc(3,3,1);
    H(1,1) = 10;
elseif P==3 | P==9
    H = sparse(diag([ 2 2 2 1/(x(4)+1)^2 2 2 2]));
elseif P==4
    H = zeros(2,2);
elseif P>=5 & P<=8
    p = Prob.user.P;
    N = Prob.user.N;
    H = spalloc( (2+N)*p, (2+N)*p, 0);
elseif P==4
    H = zeros(2,2);
elseif P==11
    H  = [ 0 0 7 6;0 0 0 0;7 0 0 0;6 0 0 0];
elseif P==12
    H  = diag([ 2; 2; 0; 0 ; 0]);
elseif P==13
    H = spalloc(11,11,0);
    H(9,10)  = -x(11);
    H(9,11)  = -x(10);
    H(10,9)  = -x(11);
    H(10,11) = -x(9);
    H(11,9)  = -x(10);
    H(11,10) = -x(9);
    H = 0.5*(H+H');
elseif P==14
    H = spalloc(5,5,0);
end

% MODIFICATION LOG
%
% 021216 ango Wrote file
% 021227 hkh  Adding 9, 10, same as 3,1; Linear constraints as nonlinear
% 030208 hkh  Adding problem 11,12
% 041020 med  Added 13, 14
% 080603 med  Switched to minlpAssign, cleaned