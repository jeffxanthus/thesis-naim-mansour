% exp_eq.m
%
% function y = exp_eq(z, Prob)
%
% Prob sends the variables: R0, R1, R2, R3, R4, P, Q, D, sig
%
% The function computes values of the equation
%
%     0 = N^4*R4/(2D^3) + N^3*R3/(2D^2) + N^2*R2/D^2 + N*R1/D) + R0, 
%
% where N=P+sqrt(Q) and R0,...,D are polynomials in z of degree less than four.
%
% Input variables:
% z        : A point or a vector.
% R0,...,D : Matrices with polynomial coefficients
% sig      : A sign which controls whether to use N = P + or - sqrt(Q).
% y        : The function value y=f(z).

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1997-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written May 26, 1997.      Last modified Nov 26, 2000.

function y = exp_eq(z, Prob)

R0  = Prob.R0;
R1  = Prob.R1;
R2  = Prob.R2;
R3  = Prob.R3;
R4  = Prob.R4;
P   = Prob.P;
Q   = Prob.Q;
D   = Prob.D;
sig = Prob.sig;

Pz = polyval(P, z);
Qz = sqrt(max(0, polyval(Q, z)));
Dz = polyval(D, z);
Dz=Dz(:);
if sig(1)
   Nz=Pz-Qz;
else
   Nz=Pz+Qz;
end
Nz=Nz(:);
z=z(:);
Rz=[polyval(R4, z), polyval(R3, z), polyval(R2, z),...
    polyval(R1, z), polyval(R0, z)];
% The row above may get problem with dimensions
Qot=Nz./Dz;

%Use nested multiplication when evaluating this polynomial-like expression.
y = (((Rz(:, 1).*Qot+Rz(:, 2)).*Nz+Rz(:, 3)).*Qot+Rz(:, 4)).*Qot+Rz(:, 5);

if length(sig)>1
   if sig(2)==1
      y=abs(y);
   elseif sig(2)==2
      if y==0
         y=-100;
      else
         y=log10(abs(y));
      end
   end
end