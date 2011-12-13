% expInit.m
%
% Find initial values for the fitting of positive sums of exponentials 
% (nonlinear least squares fitting or separable nonlinear least squares)
%
% function [Prob] = expInit (Prob, ask);
%
% INPUT:  
%   Prob   Problem structure
%   ask    1:  ask questions; 0: use defaults; -1: use values in uP
%          11: ask questions in the GUI window
%              If isempty(ask) then [if isempty(uP),ask=1, else ask=-1];
%
% OUTPUT: 
%   Prob   Problem structure
%       where expInit is setting all variables in the field ExpFit:
   
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written May 29, 1999.     Last modified Jan 16, 2003.

function [Prob] = expInit(Prob, ask)

if nargin < 2
   ask=[];
end

if isempty(ask) 
   if isempty(Prob.ExpFit) 
      ask=1; 
   else
      ask=-1; 
   end
end

% ======================================================================
% SPECIAL PART I FOR EXPONENTIAL FITTING PROBLEM

Zt=[];

y=Prob.LS.y;
t=Prob.LS.t;

if ask > 0 
   [SepAlg, p, wType, eType] = expGet1(Prob); % To get correct value for eType


   if ask==1
      fprintf('\n\n')
      fprintf('Fitting exponentials sums to data. Problem: %s\n\n',Prob.Name)
      SepAlg=-1+menu('How to solve exponential fitting problem:'...
         ,'0.  Ordinary  non linear least squares'...
         ,'1.  Separable non linear least squares. alpha solved by NNLS');

      if SepAlg
         wType=-1+menu('How to weight the residuals r(lambda,alpha):'...
         ,'0.  No weighting'...
         ,'1.  y-weighting. Weight with data y: r = (E(lambda)*alpha-y) ./ y');
      else
         wType=-1+menu('How to weight the residuals r(lambda,alpha):'...
         ,'0.  No weighting'...
         ,'1.  y-weighting. Weight with data y: r = (E(lambda)*alpha-y) ./ y'...
         ,'2.  ML-weighting.  r = (E(lambda)*alpha-y) ./ E(lambda)'...
         ,'3.  max(y,ML)-weighting.  r = (E*alpha-y) ./ max(E*alpha,y)');
      end

      ss=sprintf('Number of terms in exponential sum <= %4.0f\n\n',length(t)/2);
      p=inputSet('Give number of exponential terms: ',...
                  [1:min(8,floor(length(t)/2))], ask,ss);

      Prob=expSet1(Prob, SepAlg, p, wType, eType);
   end

elseif ask==0
   % [SepAlg, p, wType, eType] = expGet1([]);
   [SepAlg, p, wType, eType] = expGet1(Prob);
   if SepAlg
      % No ML weighting allowed if separable formulation 
      wType=min(1,wType);
      Prob.ExpFit.wType=wType;
   end
else
   [SepAlg, p, wType, eType] = expGet1(Prob);
   if SepAlg
      % No ML weighting allowed if separable formulation 
      wType=min(1,wType);
      Prob.ExpFit.wType=wType;
   end
end

N=length(t);

alpha_max=1.2*max(y);

% Test if t is equidistant

m=length(t);

if m > 2                   % Test if t is equidistant
   h=t(2:m)-t(1:m-1);
   EQt=all(abs(h-h(1)) < 1e-10);
else
   EQt=1;
end
Prob.ExpFit.t_eqDist=EQt;


[lambda, alpha, beta] = expGetLa(Prob);

if length(lambda) == p
   % OK, we have lambda
elseif ask<=0
      
   if isempty(lambda) | length(lambda) ~= p
      [alpha, lambda, beta, SQ] = exp_q(Prob);

      % Sort nonzero values in increasing orders, and zeros last.
      ip=find(lambda > 0);
      [lam,il]=sort(lambda(ip));
      il=[ip(il);find(lambda==0)];
      lambda=lambda(il);
      alpha=alpha(il);
      if ~isempty(beta)
         beta=beta(il);
      end
   end
else   % ask > 0

   [min_y j]=min(y);
   l_max=-log(min_y/alpha_max)/t(j);
   %%fprintf('\nLambda_max estimated%12.4f with alpha%12.4f\n',l_max,alpha_max);
   ss=sprintf('Lambda_max estimated: %5.4f with alpha %5.4f',l_max,alpha_max); 
   l_max=-log(2*min_y/alpha_max)/t(j);
   %%fprintf('\nLambda_max estimated%12.4f with alpha%12.4f\n\n',...
   %%          l_max,alpha_max/2);
   ss=str2mat(ss,sprintf('Lambda_max estimated: %5.4f with alpha %5.4f', ...
                           l_max,alpha_max/2)); 

   [alpha, lambda, beta, SQ ] = exp_q(Prob);

   % Sort nonzero values in increasing orders, and zeros last.
   ip=find(lambda > 0);
   [lam,il]=sort(lambda(ip));
   il=[ip(il);find(lambda==0)];
   lambda=lambda(il);
   alpha=alpha(il);
   if ~isempty(beta)
      beta=beta(il);
   end
   Qlambda=lambda(:);
   Qalpha=alpha(:);

   %%xprint(lambda,'lambda_0:'); 
   %%xprint(alpha,'alpha_0: ');
   ssTemp = 'Alpha default values:'; 
   for i=1:length(alpha)
      ssTemp = [ssTemp sprintf('  %7.6f',alpha(i))];
   end
   ss=str2mat(ss,ssTemp); 

   if ~isempty(beta)
      %%xprint(beta,'beta_0:  '); 
      ss=str2mat(ss,sprintf('beta_0: %7.6f',beta));
      %%fprintf('\n f(lambda_0,alpha_0,beta_0) = %15.8e\n',SQ);
      ss=str2mat(ss,sprintf('f(lambda_0,alpha_0,beta_0) = %9.8e',SQ)); 
   else
      %%fprintf('\n f(lambda_0,alpha_0) = %15.8e\n',SQ);
      ss=str2mat(ss,sprintf('f(lambda_0,alpha_0) = %9.8e',SQ)); 
   end

   %if 0
   % % Compute and plot transformed residual
   % ll=[0:0.05:10]';
   % El=exp(-t*ll');
   % eta=y'*El;
   % figure(2);
   % clf;
   % title('Eta plot (transformed residual in lambda [0,10] step 0.05');hold on;
   % plot(ll,eta)
   % hold on;
   % disp ('*** Transformed residual plot *** ')
   %end

   ss=str2mat(ss,sprintf('0 < Lambda(i) < Lambda(i+1)'));
   ss=str2mat(ss,sprintf('Give increasing sequence of starting lambdas.'));
   ss=str2mat(ss,sprintf(...
      'A negative value gives default initial values for all lambda;'));
   ss=str2mat(ss,sprintf(...
      'the lambda values obtained by Generalized Interpolation'));
   ss=str2mat(ss,' ');
   ssTemp = sprintf('( %7.6f )',lambda(1)); 
   lambda(1)=inputR(['Give lowest lambda starting value ' ssTemp ' : '],...
                    -1E50,30,lambda(1),ask,ss);
   if lambda(1) < 0   % Use default values, GP new
      lambda=Qlambda;
      alpha=Qalpha;   % Will be overwritten. Shall we use this instead?
   else
      for i=2:p
          ss=[];
          ss=str2mat(ss,sprintf('-----------------------------'));
          ss=str2mat(ss,' ');
          ssQuest = [sprintf('Give lambda #%d starting value ',i) ...
                     sprintf('( %7.6f ) : ',lambda(i))]; 
          lambda(i)=inputR(ssQuest,lambda(i-1),30,lambda(i),ask,ss);
      end
   end
end 

lambda=lambda(:);

if isempty(alpha) | ~length(alpha)==p
   % Estimate starting value of alpha using NNLS, alpha >=0 always
   Zt=exp(t*(-lambda)');
   %alpha=nnls1(Zt,y);
   alpha=Tnnls(Zt,y);
   if max(alpha) > 5*max(y)
      xprinte(y,'y:');
      xprinte(alpha,'alpha:');
      PrintMatrix(Zt,'Zt');
      fprintf('Max(alpha) > 5*max(y)\n')
      fprintf('Condition number in Zt*alpha = y %25.0f\n',cond(Zt))
   end
   alpha=min(2*max(y),alpha);  % Must safeguard for ill conditioning
end

if ask > 0
   if ask==1
      disp('lambda');
      disp(lambda);
      disp('alpha');
      disp(alpha);
   end
   if isempty(Zt)
      Zt=exp(t*(-lambda)');
   end
   if ask==1 % No plot when running GUI
      r=Zt*alpha-y;
      fprintf('\n *** Plot of residual at start ***\n\n')
      figure(3);
      clf;
      plot(r); title('Residual at start'); hold on;
    end
end

if SepAlg
   Prob.x_0 = lambda;
   Prob.x_L=zeros(p,1);
   Prob.x_U=[30*ones(p,1)];
   Prob.x_max=[5*ones(p,1)];
   Prob.N = p;
elseif eType < 4
   Prob.x_0 = [lambda;max(0,alpha)];
   Prob.x_L=zeros(2*p,1);
   Prob.x_U=[50*ones(p,1);5*alpha_max*ones(p,1)];
   Prob.x_max=[5*ones(p,1);5*ones(p,1)];
   Prob.N = 2*p;
else
   Prob.x_0 = [lambda;max(0,alpha);max(0,beta)];
   Prob.x_L=zeros(3*p,1);
   Prob.x_U=[30*ones(p,1);5*alpha_max*ones(2*p,1)];
   Prob.x_max=[5*ones(p,1);5*ones(2*p,1)];
   Prob.N = 3*p;
end

% Set starting values of lambda and alpha (and beta)

Prob=expSetLa(Prob, lambda, alpha, beta);

Prob.x_min=Prob.x_L;

% MODIFICATION LOG
%
% 981022  hkh  Changes for new TOMLAB design (t,y) from Prob.NLLS.
%              Delete variable flopEXPq
% 981024  hkh  If SepAlg, only wType 0 or 1 is allowed
% 981028  hkh  Prob must be used to get eType, even is ask is 0.
% 981102  hkh  Must rearrange beta if lambda not sorted correct
% 981105  hkh  Delete the assignment of Prob.f_0
% 981108  hkh  Remove all globals
% 981118  hkh  Change upper bounds to 30 for lambda
% 981208  hkh  Fix input to conform with GUI
% 981206  edr  Augment text string ss
% 990204  hkh  Must check if correct length(lambda)==p
% 020422  hkh  It is important to keep the orignal length of x_0 in Prob.N
% 030116  hkh  Change name to Tnnls
