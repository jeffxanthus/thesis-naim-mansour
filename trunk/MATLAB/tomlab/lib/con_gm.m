% con_gm.m:
%
% con_gm computes the gradient of the merit function theta(x_k)
%
% function [gtheta,OUT] = con_gm(x_k, OUT, Prob, varargin)
%
% INPUT:
% x_k     Current iterate
% OUT.c_k is the current values of the nonlinear constraints at x_k
% OUT.Ax  is the current values of the linear    constraints at x_k
%
% Merit == 0 ==> Augmented Lagrangian
% Merit == 1 ==> Non-differentiable L_1-penalty function
%
% OUTPUT:
%
% gtheta    Gradient of merit function
% OUT.g_k   Gradient vector at x_k
% OUT.cJac  Constraint gradient matrix at x_k

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1995-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written June 20, 1995. Last modified Aug 13, 2009.

function [gtheta,OUT] = con_gm(x_k, OUT, Prob, varargin)

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


if fceval >= 0    % If fceval < 0, then g_k and dc_k is computed.
   % Pick up parameters from Prob structure and cval grad/constr grad 
   Prob.g_k = []; % 070712 frhe Added to avoid fdng return Prob.g_k.
   g_k  = nlp_g( x_k(1:n), Prob, varargin{:});
   if Prob.mNonLin > 0
      dc_k = nlp_dc(x_k(1:n), Prob, varargin{:});
   else
      dc_k = [];
   end
   Ax   = OUT.Ax;
   c_k  = OUT.c_k;
else
   g_k  = Prob.g_k;
   dc_k = Prob.cJac;
   Ax   = Prob.MERIT.Ax;
   c_k  = Prob.c_k;
end
if nargout > 1
   OUT.g_k = g_k;
   OUT.cJac= dc_k;
end

if isempty(c_k) && ~isempty(dc_k)
   % If c_k not sent, recompute it
   c_k = nlp_c(x_k(1:n), Prob, varargin{:});
   c_k = c_k(:);
end
if isempty(Ax) && ~isempty(Prob.A)
   % If Ax not sent, recompute it
   Ax = Prob.A*x_k(1:n);
end

z = [x_k(1:n);Ax;c_k];

% Gradient is  +1/-1 * [I; A; dc_k]
z_k1 = z-bl;
z_k2 = bu-z;
A    = Prob.A;
mA   = size(A,1);

if merit==0 % Augmented Lagrangian
   v_k = x_k(n+1:n+m);

   gtheta = g_k;
   g_v    = zeros(m,1);

   for j=1:length(E)
       i = E(j);
       if i <= n         % Simple bounds
          gtheta(i) = gtheta(i)-(v_k(i)-pen(i)*z_k1(i));
       elseif i <= n+mA; % Linear constraints
          gtheta = gtheta-(v_k(i)-pen(i)*z_k1(i))*A(i-n,:)';
       else              % General constraints
          k = i-n-mA;
          gtheta = gtheta-(v_k(i)-pen(i)*z_k1(i))*dc_k(k,:)';
       end
       g_v(i) = -v_k(i)/pen(i);
   end
   for j = 1:length(I)
       i = I(j);
       if min(z_k1(i),z_k2(i))*pen(i) <= v_k(i)
          if i <= n
             if z_k1(i) <= z_k2(i)
                gtheta(i) = gtheta(i)-(v_k(i)-pen(i)*z_k1(i));
                g_v(i) = -z_k1(i);
             else
                gtheta(i) = gtheta(i)+(v_k(i)+pen(i)*z_k2(i));
                g_v(i) = -z_k2(i);
             end
          elseif i <= n+mA;
             if z_k1(i) <= z_k2(i)
                gtheta = gtheta-(v_k(i)-pen(i)*z_k1(i))*A(i-n,:)';
                g_v(i) = -z_k1(i);
             else
                gtheta = gtheta-(v_k(i)+pen(i)*z_k2(i))*A(i-n,:)';
                g_v(i) = -z_k2(i);
             end
          else
             k = i-n-mA;
             if z_k1(i) <= z_k2(i)
                gtheta = gtheta-(v_k(i)-pen(i)*z_k1(i))*dc_k(k,:)';
                g_v(i) = -z_k1(i);
             else
                gtheta = gtheta-(v_k(i)+pen(i)*z_k2(i))*dc_k(k,:)';
                g_v(i) = -z_k2(i);
             end
          end
       else
          g_v(i) = -v_k(i)/pen(i);
       end
   end
   gtheta = [gtheta;g_v];
elseif merit == 1 
   gtheta = g_k;

   for j=1:length(E)
       i = E(j);
       if i <= n
          if z_k1(i) >= 0
             gtheta(i) = gtheta(i)+pen(i);
          else
             gtheta(i) = gtheta(i)-pen(i);
          end
       elseif i <= n+mA;
          if z_k1(i) >= 0
             gtheta = gtheta+pen(i)*A(i-n,:)';
          else
             gtheta = gtheta-pen(i)*A(i-n,:)';
          end
       else
          k = i-n-mA;
          if z_k1(k) >= 0
             gtheta = gtheta+pen(i)*dc_k(k,:)';
          else
             gtheta = gtheta-pen(i)*dc_k(k,:)';
          end
       end
   end
   for j=1:length(I)
       i = I(j);
       if min(z_k1(i),z_k2(i)) <= 0
          if i <= n
             if z_k1(i) <= z_k2(i)
                gtheta(i) = gtheta(i)-pen(i);
             else
                gtheta(i) = gtheta(i)+pen(i);
             end
          elseif i <= n+mA;
             if z_k1(i) <= z_k2(i)
                gtheta = gtheta-pen(i)*A(i-n,:)';
             else
                gtheta = gtheta+pen(i)*A(i-n,:)';
             end
          else
             k=i-n-mA;
             if z_k1(k) <= z_k2(k)
                gtheta = gtheta-pen(i)*dc_k(k,:)';
             else
                gtheta = gtheta+pen(i)*dc_k(k,:)';
             end
          end
       end
   end
end

% MODIFICATION LOG:
%
% 981014  hkh  Bound bug k should be i in test: if z_k1(i) <= z_k2(i)
% 981016  hkh  Index (i) was missing for single variable constraints in one
%              case: gtheta(i)=gtheta(i)+(v_k(i)+pen(i)*z_k2(i));
% 981017  hkh  Ax was wrong. Must be recomputed like c_k and dc_k.
%              Changed input/output using OUT
% 981110  hkh  Change field dc_k to cJac
% 001218  hkh  Change to use transposed format for dc_k, rows for constraints
% 040506  hkh  Safeguard c_k, changing to a column vector
% 070712  frhe Reset Prob.g_k if recomputation.
% 080607  hkh  Avoid call to nlp_dc if no nonlinear constraints
% 090813  med  mlint check