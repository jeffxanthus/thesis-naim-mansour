% con_fm.m:
%
% con_fm computes the merit function theta(x_k) for test problem P
%
% function [theta, OUT] = con_fm(x_k, OUT, Prob, varargin)
%
% Merit == 0 ==> Augmented Lagrangian
% Merit == 1 ==> Non-differentiable L_1-penalty function
%
% OUT.f_k   Function value  at x_k
% OUT.c_k   Constraint vector at x_k
% OUT.Ax    A*x_k  (A = linear constraints)

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.1.0$
% Written June 20, 1995. Last modified June 7, 2004.

function [theta, OUT] = con_fm(x_k, OUT, Prob, varargin)

% Pick up variables from Prob struct
merit = Prob.MERIT.merit;
pen   = Prob.MERIT.pen;
fceval= Prob.MERIT.fceval;
bl    = Prob.MERIT.bl;
bu    = Prob.MERIT.bu;
E     = Prob.MERIT.E;
I     = Prob.MERIT.I;

m = length(bl);
n = max(length(Prob.g_k),size(Prob.A,2));

if fceval >= 0    % If fceval < 0, then f_k and c_k is computed.
   % Pick up parameters from Prob structure and eval func/constr 
   f_k = nlp_f(x_k(1:n), Prob, varargin{:} );
   if Prob.mNonLin > 0
      c_k = nlp_c(x_k(1:n), Prob, varargin{:} );
      c_k = c_k(:);
   else
      c_k = [];
   end
   if isempty(Prob.A)
      Ax=zeros(0,1);
   else
      Ax=Prob.A*x_k(1:n);
   end
else
   f_k=Prob.f_k;
   c_k=Prob.c_k;
   Ax=Prob.MERIT.Ax;
end

if nargout > 1
   OUT.f_k = f_k;
   OUT.c_k = c_k;
   OUT.Ax  = Ax;
end

z      = [x_k(1:n);Ax;c_k];
z_k    = zeros(m,1);
z_k(I) = min(z(I)-bl(I),bu(I)-z(I));
z_k(E) = z(E)-bl(E);

if merit==0 % Augmented Lagrangian
   v_k = x_k(n+1:n+m);

   if isempty(I)
      k = [];
   else
      k = z_k(I).*pen(I) <= v_k(I);
   end
   J = [E;I(find(k))];
   K = I(find(~k));

   theta = f_k-sum(z_k(J).*(v_k(J)-0.5*pen(J).*z_k(J)))-...
           0.5*sum(v_k(K).^2./pen(K));

elseif merit==1 % Non-differentiable L_1-penalty function

   theta = f_k+sum(pen(E).*abs(z_k(E)))-...
               sum(pen(I).*min(z_k(I),0));
end
%if 0
%   fprintf('Merit function %15.6f %15.6f %25.15f\n',x_k(1:2),theta);
%end

% MODIFICATION LOG:
%
% 981017  hkh  Ax was wrong. Must be recomputed like c_k. 
%              Ax sent back in OUT like c_k
% 990627  hkh  Add OUT as input
% 040506  hkh  Safeguard c_k, changing to a column vector
% 080607  hkh  Avoid call to nlp_c if no nonlinear constraints
