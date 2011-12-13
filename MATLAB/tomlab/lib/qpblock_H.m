% function H = qpblock_H(x, Prob)
%
% Compute Hessian for quadratic factorized problem
%
% x      Point x where H(x) is evaluated
% Prob   Problem structure
% H      Hessian matrix, H(x) = F + Fb.out' * Fb.inn * Fb.out (for i=1...p).
%        or H(x) = F + Fb.out' * Fb.out (for i=1...p). See qpblock_f.

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2006-2009 by Tomlab Optimization Inc., Sweden. $Release: 7.4.0$
% Written Aug 21, 2006.   Last modified Dec 22, 2006.

function H = qpblock_H(x, Prob)
n=length(x);
if isempty(Prob.QP.F)
    H = sparse(n,n);
else
    if Prob.QP.case == 2 | Prob.QP.case == 4
        H = spdiags(Prob.QP.F,0,n,n);
    else 
        H = Prob.QP.F;
    end
end
if Prob.QP.case == 1 | Prob.QP.case == 3
    for i=1:Prob.QP.mQP
        H = H + Prob.QP.Fb(i).out'*Prob.QP.Fb(i).inn*Prob.QP.Fb(i).out;
    end
    H = 0.5*(H+H');
else
    for i=1:Prob.QP.mQP
        H = H + Prob.QP.Fb(i).out'*Prob.QP.Fb(i).inn*Prob.QP.Fb(i).out;
    end
    H = 0.5*(H+H');
end

% MODIFICATION LOG
%
% 060821  med  Written
% 060821  hkh  Removed parenthesis in 3-matrix computations
% 091223  hkh  Must safe guard and symmetrize H.
