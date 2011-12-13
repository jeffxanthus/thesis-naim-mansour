% Compute relative or absolute constraint violations, return L1-norm
% Optionally return all violations in one vector
%
% function [h_L1, h] = consviolation(Result,PriLev,absviol)
%
% INPUT:
%
% Result  Output result structure from Tomlab
% PriLev  =1, Print h_L1 and if any violations in lower and upper bounds of
%         variables, print number of violations of these
% absviol if == 1, compute absolute violations
%         Default absviol == 0, compute relative violations
%
% OUTPUT:
% h_L1     L1 of violations, i.e. sum of abs of all (relative) violations
% h        (Relative) violations in variables x, linear constraints A*x
%          and nonlinear constraints c(x)

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Dec 4, 2003.   Last modified Aug 13, 2009.

function [h_L1, h] = consviolation(Result,PriLev,absviol)

if nargin < 3
   absviol = 0;
   if nargin < 2
      PriLev = [];
      if nargin < 1
         fprintf('consviolation: ERROR! No input parameters\n');
         error('consviolation(Result,PriLev) is the call');
      end
   end
end
if isempty(PriLev), PriLev = 1; end

h_L1 = NaN;
h    = [];

if isfield(Result,'Prob')
   Prob = Result.Prob;
else
   if isempty(Result)
      fprintf('\nconsviolation: Empty Result. \n\n');
      fprintf('Can not compute constraint violation\n');
   else
      fprintf('\n\n\nconsviolation: ');
      fprintf('Result.Prob not defined\n\n\n');
      fprintf('Can not compute constraint violation\n');
   end
   return
end

x_k = Result.x_k;

if ~isempty(x_k) 
   if size(Prob.A,1) > 0 
      Ax=DefPar(Result,'Ax',[]);
      if isempty(Ax)
         Ax=Prob.A*x_k(:,1); 
      end
   else
      Ax=zeros(0,1);
   end

   if isempty(Result.c_k)
      zz=[x_k(:,1);Ax]; 
      bl=[Prob.x_L;Prob.b_L];
      bu=[Prob.x_U;Prob.b_U];
   else
      zz=[x_k(:,1);Ax;Result.c_k(:,1)]; 
      mL=length(Result.c_k);
      if isempty(Prob.c_L)
         bl=[Prob.x_L;Prob.b_L;-Inf*ones(mL,1)];
      else
         bl=[Prob.x_L;Prob.b_L;Prob.c_L];
      end
      if isempty(Prob.c_U)
         bu=[Prob.x_U;Prob.b_U;Inf*ones(mL,1)];
      else
         bu=[Prob.x_U;Prob.b_U;Prob.c_U];
      end
   end
   zz = real(zz);
   if length(bl)==length(zz) & length(bu)==length(zz)
       if nargout > 1
           if absviol == 1
               h =-min(0,zz-bl)-min(0,bu-zz);
           else
               h =-min(0,(zz-bl)./max(1,abs(bl)))-min(0,(bu-zz)./max(1,abs(bu)));
           end
           h_L1 = sum(h);
       else
           if absviol == 1
               h_L1 = norm(min(0,zz-bl),1) + norm(min(0,bu-zz),1);
           else
               h_L1 = norm(min(0,(zz-bl)./max(1,abs(bl))),1) + ...
                   norm(min(0,(bu-zz)./max(1,abs(bu))),1);
           end
       end
   else
       fprintf('Conflicting dimensions on bounds and constraint values\n');
       fprintf('Length(bl) %d Length(bu) %d Length([x,Ax,c_k]) %d\n',...
           length(bl),length(bu),length(zz));
       NN=min([length(bl),length(bu),length(zz)]);
       if nargout > 1
           if absviol == 1
               h= max(0,bl(1:NN)-zz(1:NN))+max(0,zz(1:NN)-bu(1:NN));
           else
               h=  max(0,(bl(1:NN)-zz(1:NN))./max(1,abs(bl(1:NN)))) + ...
                   max(0,(zz(1:NN)-bu(1:NN))./max(1,abs(bu(1:NN))));
           end
           h_L1 = sum(h);
       else
           if absviol == 1
               h_L1 = norm(max(0,bl(1:NN)-zz(1:NN)),1) + ...
                   norm(max(0,zz(1:NN)-bu(1:NN)),1);
           else
               h_L1=norm(max(0,(bl(1:NN)-zz(1:NN))./max(1,abs(bl(1:NN)))),1) +...
                   norm(max(0,(zz(1:NN)-bu(1:NN))./max(1,abs(bu(1:NN)))),1);
           end
       end
   end
   if PriLev >0 & ~isempty(h_L1) % & ~(h_L1==0)
      fprintf('                      sum(|constr|)%26.18f',h_L1);
      fprintf('\n');
   end

   if PriLev > 0
      xTol       = Prob.optParam.xTol;
      vL=sum(Prob.x_L > x_k(:,1)+xTol);
      vU=sum(Prob.x_U < x_k(:,1)-xTol);
      if vL+vU > 0
         fprintf(' ==>  Number of variables violating lower bound%4d. ',vL);
         fprintf(' Number of variables violating upper bound%4d',vU);
         fprintf('\n');
      end
   end
else
   fprintf('No optimal solution Result.x_k found\n');
   h_L1 = NaN;
   h    = [];
end

% MODIFICATION LOG:
%
% 031204  hkh  Written
% 040324  hkh  typo bug
% 080417  hkh  Change to relative violations, keep option to compute absolute
% 090813  med  mlint check