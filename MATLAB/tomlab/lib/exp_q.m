% exp_q.m
%
% function [a, b, c, SQ ] = exp_q(Prob)
%
% The function tries to solve the problem of fitting different
% kinds of exponential sums (ES) to data (t, y).
%
% (ES1)    a(1)*exp(-b(1)*t)+...+a(p)*exp(-b(p)*t)
% (ES2)    a(1)(1-*exp(-b(1)*t))+...+a(p)*(1-exp(-b(p)*t))
% (ES3)    a(1)*t*exp(-b(1)*t)+...+a(p)*t*exp(-b(p)*t)
% (ES4)    (a(1)*t+c(1))*exp(-b(1)*t)+...+(a(p)*t+c(p))*exp(-b(p)*t)
%
% The program exp_geo is developed for data equidistant in t. So first
% data must be transformed to (t, y) equidistant in t.
%
% The program is divided into some natural parts. 
% 1. Preprocess data series for the relevant kind of ES so that the
%    method for making data equidistant is applicable.
% 2. Make data equidistant in t.
% 3. Compute different S and loop over this to find estimates
%    for nonlinear variables b (by z from exp_geo). 
% 4. Compute linear variables a by solving a NNLS problem, using Tnnls.
% 5. Make a choice between suggested solution candidates (a, b) 
%    using a weighted LS-criterion.
%
% INPUT:  Problem structure Prob.
%
%  In Prob.LS fields used are:
%  t, y  : a time series (t, y). 
%
%  In Prob.ExpFit use fields:
%  p      : length(lambda) asked for.
%  eType:   Choice of algorithm depending on type of exponential sum ES 
%  t_eqd:   True if the series is equidistant 
%  wType    Type of weighting 
%  dType    Differentiation formula for ES1. 
%  geoTyp   Type of equation, used by exp_geo. 
%  qType    Length q of partial sums. 
%  sigType  Sign to use in (P +/- sqrt(Q))/D in exp_geo for p==3 or 4.
%
% Output variables:
%  a, b : Linear and nonlinear coefficients in the exponential sum. 
%  SQ   : Sum of squares for the estimate.
%
% Algorithms:
%
% Choice of part of a series:
%   The main rule is to choose a series so that |y(j)| is big. 
%   The reason for this is to maximize the ratio of signal to noise.
%
% Choice of length of a part of a series:
%   Partials series of lengths 2*p*q with q=1:4 are chosen if possible 
%   (the number of data points must be at least 2*p*q). 
%   If the noise is low, q=1 could give good results.
%   And for decaying series q=4 will make the series so long that fast 
%   terms seem to vanish too quickly.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Aug 17, 1999.  Last modified Aug 14, 2006.

function [a, b, c, SQ] = exp_q(Prob)

t=Prob.LS.t;
y=Prob.LS.y;

SSQ=[];         % Sum of squares. Used as a criterion of choice.
A=[];B=[];C=[]; % Storage of variables.
n0=length(t);

%Initialization of parameters
p=Prob.ExpFit.p;
eType=Prob.ExpFit.eType;
equid=Prob.ExpFit.t_eqDist; %Equidistant if equid==1
wType=Prob.ExpFit.wType; %Weight
dType=Prob.ExpFit.dType; %Differentiation scheme
x0Type=Prob.ExpFit.x0Type; %Type of starting value

XgeoType=Prob.ExpFit.geoType; % geoType
XqType  =Prob.ExpFit.qType;   % qType 
XsigType=Prob.ExpFit.sigType; % sigType 

sumType=Prob.ExpFit.sumType;
if sumType==0
   sType=1:3;
else
   sType=sumType;
end

%Some special features for different eType.
eType=max(1,eType); %if eType<1.
if eType==1      %Some special features for eType==1
   dType=0;
   if XgeoType==0; geoType=[0, 1]; %Both heuristic and analytic algorithm.
   elseif XgeoType==1; geoType=1; %Analytic algorithm.
   else geoType=0;  %Heuristic algorithm.
   end   
   if p(1)~=4;
      geoType=0;  %exp_geo uses a combination of 
      %analytic and heuristic algorithm.
   end
   ESt=[n0;t]; ESy=[n0;y];
elseif eType==2       % Some special features for eType==2
   if equid==0;       %If not equidistant,
      if dType==0;    %and no differentiantion is asked for,
         dType=[1:3]; %then choose some differentiation.
      end      
   else
      dType=0;  %use no differentiation
   end
   if dType==0 & (x0Type==-4 | x0Type==-5)
      dType=1:3; %Differentiated Prony. 
      %else if x0Type==-3 non differentiated Prony.
   end   
   %Column 1 is original data and   
   %column 2, 3 and 4 in ESt and ESy is forward, backward and
   %central difference quotients.
   ESt=[ [n0;t], [n0-1;t(1:n0-1);0], [n0-1;t(2:n0);0],...
         [n0-2;t(2:n0-1); 0;0] ];
   diff1=(y(2:n0)-y(1:n0-1))./(t(2:n0)-t(1:n0-1));
   diff2=(y(3:n0)-y(1:n0-2))./(t(3:n0)-t(1:n0-2));
   ESy=[ [n0;y], [n0-1;diff1(:);0], [n0-1;diff1(:);0],...
         [n0-2;diff2(:); 0;0] ];
   Mt=ESt; My=ESy;
   geoType=0;  %Analytic algorithm for p==3 and heuristic when p==4.
   %sigType=3;
elseif eType==3      % Some special features for eType==3
   geoType=0;  
   dType=0;
   ESt=[n0;t]; ESy=[n0;y];
elseif eType==4      % Some special features for eType==4
   geoType=0;  
   dType=0;
   ESt=[n0;t]; ESy=[n0;y];
end

Mt=[]; My=[];
if equid==0
   if eType==1 | eType==2;
      % Construction of equidistant times.
      for i1=1:length(dType);      
         n1=ESt(1,dType(i1)+1);
         tprel=ESt(2:n1+1, dType(i1)+1);
         Yprel=ESy(2:n1+1, dType(i1)+1);
         delta=tprel(2)-tprel(1);
         te=[tprel(1):delta:tprel(n1)]';
         n2=length(te);
         %ML-reconstruction of interpolating equidistant data Ye.
         Ye=Yprel(1:2);
         for i2=2:n1-1;
            if length(Ye)<n2;
               %This if-construction prevents out of index in te if
               %there is an interval [tprel(k:k+1)] beyond [te(n2-1:n2)].
               tsubint=[te(length(Ye)+1):delta:tprel(i2+1)]';
               if length(tsubint)>0
                  %This if-construction is needed if next te(j) is 
                  %outside tprel(i2:i2+1).
                  if Yprel(i2)*Yprel(i2+1)<0;
                     %This if-construction is needed as the algorithm 
                     %does not cope a change of sign in Yprel. In such
                     %case use linear interpolation.
                     Ynew=( (Yprel(i2+1)-Yprel(i2))/(tprel(i2+1)-tprel(i2)) ) *...
                        ( tsubint-tprel(i2) ) + Yprel(i2);
                  else
                     be=log( Yprel(i2+1)/Yprel(i2) )/( tprel(i2)-tprel(i2+1) );
                     ae=Yprel(i2)*exp(be*tprel(i2));
                     Ynew=ae*exp(-be*tsubint);
                  end
                  Ye=[Ye;Ynew(:)];
               end
            end   
         end
         Mt(1,dType(i1)+1)=n2;       My(1,dType(i1)+1)=n2;
         Mt(2:n2+1,dType(i1)+1)=te;  My(2:n2+1,dType(i1)+1)=Ye;
      end %for ...dType   
   elseif eType==3;
      % Construction of equidistant times.
      n1=ESt(1,1);
      tprel=ESt(2:n1+1, 1);
      Yprel=ESy(2:n1+1, 1);
      delta=tprel(2)-tprel(1);
      te=[tprel(1):delta:tprel(n1)]';
      n2=length(te);
      %Construction of interpolating equidistant data Ye.
      Ye=Yprel(1:2);
      for i2=2:n1-1;
         if length(Ye)<n2;
            %This if-construction prevents out of index in te if
            tsubint=[te(length(Ye)+1):delta:tprel(i2+1)]';
            if length(tsubint)>0
               %This if-construction is needed if next te(j) is 
               %outside tprel(i2:i2+1).
               if Yprel(i2)*Yprel(i2+1)<0;
                  %This if-construction is needed as the algorithm 
                  %does not cope a change of sign in Yprel. In such
                  %case use linear interpolation.
                  Ynew=( (Yprel(i2+1)-Yprel(i2))/(tprel(i2+1)-tprel(i2)) ) *...
                     ( tsubint-tprel(i2) ) + Yprel(i2);
               else
                  %be=log( Yprel(i2+1)*tprel(i2)/tprel(i2+1)/Yprel(i2) )/...
                  %   (tprel(i2+1)-tprel(i2));
                  %ae=( Yprel(i2+1)/tprel(i2+1) )*exp( be*tprel(i2+1) ); 
                  %Ynew=ae*tsubint.*exp(-be*tsubint);
                  bas=Yprel(i2)*tprel(i2+1)/Yprel(i2+1)/tprel(i2);
                  pot=(tprel(i2)-tsubint)/(tprel(i2+1)-tprel(i2));
                  amp=tsubint*Yprel(i2)/tprel(i2);
                  Ynew=amp.*bas.^pot;
               end
               Ye=[Ye;Ynew(:)];
            end
         end   
      end
      Mt(1,1)=n2;       My(1,1)=n2;
      Mt(2:n2+1,1)=te(:);  My(2:n2+1,1)=Ye(:);
   elseif eType==4;
      % Construction of equidistant times.
      n1=ESt(1,1);
      tprel=ESt(2:n1+1, 1);
      Yprel=ESy(2:n1+1, 1);
      delta=tprel(2)-tprel(1);
      te=[tprel(1):delta:tprel(n1)]';
      n2=length(te);      
      %Construction of interpolating equidistant data Ye.
      Ye = interp1(t,y,te, 'spline');
      %The alternatives  'linear', 'cubic' and 'nearest' in interp1
      %seemed not to be as good as 'spline'.
      %Neither the alternative below works good. It uses interpolation
      %by estimating y' and solving y'/y=1/(t+c/a)-b.
      %if 0
      %   n0=length(t);
      %   ESt=[n0;t(:)]; ESy=[n0;y(:)];
      %   intn=floor(n1/2);
      %   if intn~=n1/2; %If even number of data points
      %      Yprel=y(1:n1-1); tprel=t(1:n1-1); %Omit last "odd point".
      %   else
      %      Yprel=y; tprel=t; 
      %   end
      %   Ye=Yprel(1:2);      
      %   for i1=1:intn-1;
      %      tsubint=te(sum(te<=t(2*i1)):sum(te<=t(2*i1+2)));
      %      %Compute a, b, c.
      %      e=(ESt(i1)-2*tprel(i1:i1+1)+ESt(i1+2))/...
      %         (ESt(i1)-ESt(i1+1))/(ESt(i1+1)-ESt(i1+2));
      %      f = (2*tprel(i1:i1+1)-ESt(i1)-ESt(i1+1))/...
      %         (ESt(i1)*ESt(i1+1)-ESt(i1)*ESt(i1+2)+ESt(i1+2)^2-...
      %          ESt(i1+1)*ESt(i1+2));
      %      d = (2*tprel(i1:i1+1)-ESt(i1+1)-ESt(i1+2))/...
      %         ((ESt(i1)-ESt(i1+2))*(ESt(i1)-ESt(i1+1)));
      %      yprim=d*ESy(i1)+e*ESy(i1+1)+f*ESy(i1+2);         
      %      %Choose e.g. forward quotient below.
      %      sq=(ESt(i1)-ESt(i1+1)-4/(yprim(1)/ESy(i1)-yprim(2)/ESy(i1+1)))*...
      %         (ESt(i1)-ESt(i1+1));
      %      g=-0.5*( ESt(i1)+ESt(i1+1)-sqrt(sq));
      %      b=yprim(2)/ESy(i1+1)-1/(ESt(i1+1)+g);
      %      a=( yprim(1)+b*ESy(i1) ) /exp(-b*ESt(i1));
      %      c=g*a;
      %      Ynew=(a*tsubint+c).*exp(-b*tsubint);
      %      Ye=[Ye;Ynew(:)];
      %   end
      %end  %if 0      
   else
      disp('This eType not implemented in exp_q');
   end
else 
   Mt(1,1)=n0;  My(1,1)=n0;
   Mt(2:n0+1)=t;  My(2:n0+1)=y;
   if eType==2;
      if dType==0 & (x0Type==-4 | x0Type==-5)
         dType=1:3; %Differentiated Prony. 
         %else if x0Type==-3 non differentiated Prony.
      end   
      %Column 1 is original data and   
      %column 2, 3 and 4 in ESt and ESy is forward, backward and
      %central difference quotients.
      Mt=[ [n0;t], [n0-1;t(1:n0-1);0], [n0-1;t(2:n0);0],...
            [n0-2;t(2:n0-1); 0;0] ];
      diff1=(y(2:n0)-y(1:n0-1))./(t(2:n0)-t(1:n0-1));
      diff2=(y(3:n0)-y(1:n0-2))./(t(3:n0)-t(1:n0-2));
      My=[ [n0;y], [n0-1;diff1(:);0], [n0-1;diff1(:);0],...
            [n0-2;diff2(:); 0;0] ];
   end   
end   

if (x0Type==-3 | x0Type==-4) & exist('mprony','file')
   if eType==1
      b = feval('mprony',y,t,p(1));
      b=max(real(-b(:)),0);
      Zt=exp(t*(-b)');
      %a=nnls1(Zt, y);
      a=Tnnls(Zt, y);
      c=[]; SQ=[];
      ys=Zt*a;
   elseif eType==2
      %if x0Type==-3 | x0type==-4
         if dType==0
            b = feval('mprony',y,t,p(1)+1);
            b=sort(max(0,-real(b)));
            b=b(2:length(b)); b=b(:);
         else
            [b,a] = feval('mprony',My(2:My(1,4), 4),Mt(2:Mt(1,4)),p(1));
            b=sort(max(0,-real(b(:))));
         end
      %end
      Zt=1-exp(t*(-b)');
      %a=nnls1(Zt, y);
      a=Tnnls(Zt, y);
      c=[]; SQ=[];
      ys=Zt*a;
   elseif eType>2
      [b,a] = feval('mprony',y,t,p(1));
      b=sort(max(-real(b(:)),0));
      for i=1:p(1);
         sols(i)=mean(b(2*i-1:2*i));
      end   
      b=sols(:);
      Zt=exp(t*(-b)');
      tZt=(t*ones(1,p(1))).*Zt;
      if eType==3;
         %a=nnls1(tZt, y);
         a=Tnnls(tZt, y);
         c=[];
         ys=tZt*a;
      elseif eType==4;
         %a=nnls1([tZt,Zt],y);
         a=Tnnls([tZt,Zt],y);
         c=a(p+1:p+p); a=a(1:p);
         ys=tZt*a+Zt*c;
      end
   end
   if wType==0;
      W=ones(size(y));
   elseif wType==1;
      W=abs(y);
   elseif wType==2;
      W=abs(ys);
   elseif wType==3;
      W=max(abs(ys), abs(y));
   end   
   for i6=1:length(W);
      if W(i6)==0; W(i6)=1; end 
      %Avoid division by zero by using weigth=1.
   end            
   SQ=norm((ys-y)./W, 2);
else
%Preparations for exp_geo
p(2)=eType;
if eType==2
   if dType(1)>0;
      p(2)=1;
      %Entering derivatives of series to exp_geo ->treat as eType==1.
   end
end  
p(3)=1; %=z_max=1 if positive expsums.
%Compute estimates.
A=[]; B=[]; C=[]; SSQ=[];

for i1=1:length(geoType);
   %qmax = Maximal length of partial sum. What to do if qmax==0 occurs?
   %If there are restrictions in sign of parameters (a, b), this case is
   %not irrelevant (but rare), although it means that there are fewer data
   %points then parameters.   
   peff=min(p(1), 3)+max(geoType);
   if (x0Type==-5 & eType>2); peff=2*p(1); end %peff dependent on algorithm.
   %For p==4, geoType(i)==1 will give a numerical estimate and else a heuristic estimate.
   for i2=1:length(dType);
      if dType(1)==0
         n2=n0; te=t; Ye=y;
      else
         n2=Mt(1, dType(i2)+1);
         te=Mt(2:n2+1, dType(i2)+1);
         Ye=My(2:n2+1, dType(i2)+1);
      end         
      if eType==1
         qmax=floor(n2/(2*peff));
      elseif eType==2
         qmax=floor(n2/(2*peff+1));
      elseif eType==3 | eType==4
         if x0Type==0
            qmax=floor(n2/(3*peff)); %For GInterp.
         elseif x0Type==-5
            qmax=floor(n2/(2*peff)); %For GProny.
         end
      else
         %Not implemented yet.
      end      
      if XqType==0; qType=[1: qmax]; 
      else qType=min(XqType, qmax);
      end
      for i3=1:length(qType);
         %Divide the series into partial sums S(i4). 
         q=qType(i3);
         if p(2)==1  %Can't use eType here as (te, Ye) might 
            %not be of eType.
            Slen=2*peff;
         elseif p(2)==2;            
            if x0Type==0
               Slen=2*peff+1; %For GInterp.
            elseif x0Type==-5
               Slen=2*peff; %For GProny
            end            
         elseif p(2)==3 | p(2)==4;
            if x0Type==0
               Slen=3*peff; %For GInterp.
            elseif x0Type==-5
               Slen=2*peff; %For GProny
            end            
         end
         S=zeros(Slen, 1);
         for i5=1:length(sType)
            if sType(i5)==1
               for i4=1:Slen            
                  S(i4)=sum(Ye(1+(i4-1)*q:i4*q));
               end
            elseif sType(i5)==2
               for i4=1:Slen            
                  S(i4)=sum(Ye(i4:i4+q-1));
               end
            elseif sType(i5)==3
               for i4=1:Slen            
                  S(i4)=sum(Ye(i4+([0:q-1]*2*peff)));
               end
            end
            if XsigType==0 & (p(1)==3 | peff==4) & eType==1
               sigType=1:8; 
            else
               sigType=max(XsigType, 1); 
            end
            for i4=1:length(sigType)
               p(4)=sigType(i4);
               %Compute estimates.
               %The calls exp_geo(p, te, Ye); and [z, z_lim] = exp_geo(p, S); are
               %possible. (With either data series or partial sum).
               if x0Type==0
                  [z, z_lim] = exp_geo(p, S);
                  % NOTE! Must add this entry for eType == 4 !!!
                  if isempty(z)
                     [z, z_lim] = exp_geo(p, te, Ye);
                  end
               elseif x0Type==-5
                  %z=min(max(eps,G_prony(S, eType,Prob)), p(3));
                  %z=sort(z);
                  if eType==1
                     b = feval('mprony',S,[te(1), te(2)],p(1));
                     b=max(real(-b(:)),0);
                  elseif eType==2
                     if dType==0
                        b = feval('mprony',S,[te(1), te(2)],p(1)+1);
                        b=sort(max(0,-real(b)));
                        b=b(2:length(b)); b=b(:);
                     else
                        [b] = feval('mprony',S,[te(1), te(2)],p(1));
                        b=sort(max(0,-real(b(:))));
                     end
                  elseif eType>2
                     b = feval('mprony',S,[te(1), te(2)],p(1));
                     b=sort(max(-real(b(:)),0));
                     for i=1:p(1);
                        sols(i)=mean(b(2*i-1:2*i));
                     end   
                     b=sols(:);
                  end
                  z=exp(-b*(te(2)-te(1)));
               end
               z=z(:);               
               if length(z)~=p(1);
                  disp('length(z)~=p(1)');
               end               
            
               % COMPUTE CRITERION FOR ESTIMATES.
               t_step=(te(n2)-te(1))/(n2-1); %Length of equidistant steps.
               b=-log(z)/(t_step);
               if sType(i5)==1; b=b/q; end
               Zt=exp(t*(-b)');
               if eType==1;
                  %a=nnls1(Zt, y);
                  a=Tnnls(Zt, y);
                  ys=Zt*a;
               elseif eType==2;
                  %a=nnls1(1-Zt, y);
                  a=Tnnls(1-Zt, y);
                  ys=(1-Zt)*a;
               elseif eType==3;
                  tZt=(t*ones(1,p(1))).*Zt;
                  %a=nnls1(tZt, y);
                  a=Tnnls(tZt, y);
                  ys=tZt*a;
               elseif eType==4;
                  tZt=(t*ones(1,p(1))).*Zt;
                  %a=nnls1([tZt,Zt],y);
                  a=Tnnls([tZt,Zt],y);
                  c=a(p+1:p+p); a=a(1:p);
                  ys=tZt*a+Zt*c;
                  %ys=[(t*a'+ones(size(t))*c').*exp(t*(-b)')]';  
               else
                  disp('Not ready yet for eType>4 in exp_q.m');
               end
               if wType==0;
                  W=ones(size(y));
               elseif wType==1;
                  W=abs(y);
               elseif wType==2;
                  W=abs(ys);
               elseif wType==3;
                  W=max(abs(ys), abs(y));
               end   
               for i6=1:length(W);
                  if W(i6)==0; W(i6)=1; end 
                  %Avoid division by zero by using weigth=1.
               end            
               Wr=(ys-y)./W;
               SSQ=[SSQ, Wr'*Wr/2];  B=[B, b(:)]; A=[A, a(:)];  
               if eType==4; C=[C, c(:)]; end;
            end %for i4=1:length(sigType);
         end %for i5=1:length(sType);
      end %for i3=1:length(qType);
   end %for i2=1:length(dType);
end %for i1=1:length(geoType);

%Criterion of choice is the smallest sum of squares.
[SQ, i]=min(SSQ);

a=A(:, i); b=B(:, i); 

% HKH ADDED 2002-08-02
if any(b==0)
   % Try to find some other solution with nonzero(b)
   while ~all(isinf(SSQ)) & any(b==0)
      SSQ(i) = Inf;
      [SQ, i]=min(SSQ);
      a=A(:, i); b=B(:, i); 
   end
   if any(b==0)
      % Create some kind of nonzero b  solution
      bSum=sum(b);
      if bSum == 0, bSum = 1; end
      switch p(1)
        case 1
             b(1) = 1;
        case 2
             b(1) = bSum/4;
             b(2) = 3*bSum/4;
        case 3
             b(1) = bSum/9;
             b(2) = 3*bSum/9;
             b(3) = 5*bSum/9;
        case 4
             b(1) = bSum/16;
             b(2) = 3*bSum/16;
             b(3) = 5*bSum/16;
             b(4) = 7*bSum/16;
        otherwise
             i    = find(b==0);
             b(i) = [1:length(i)]*bSum/10;
      end
   end
end
if eType==4; c=C(:, i); 
else c=[];
end
end %if x0Type==-3

% end % of routine exp_q

% MODIFICATION LOGS
%
% 980529   hkh   Made modifications for new TOMLAB using structure Prob.
% 980603   hkh   Merged conflicting versions of exp_q
% 980608   hkh   Merged conflicting versions of exp_q
% 981021   hkh   Changed access to (t,y) in Prob
% 981022   hkh   Delete variable flopEXPq
% 011205   hkh   Must safeguard call to exp_geo if eType == 4


%
% =======================================================================
%        exp_geo.m
% =======================================================================
%
% function [z, z_lim] = exp_geo(p, t, y)
%  or
% function [z, z_lim] = exp_geo(p, S)
%
% The function tries to solve the problem of fitting a function 
%
% (eType=1)    a(1)*exp(-b(1)*t)+...+a(p)*exp(-b(p)*t) 
% (eType=2)    a(1)*(1-exp(-b(1)*t))+...+a(p)*(1-exp(-b(p)*t))  
%
% to empirical data (t, y) equidistant in t.
%
% It uses expressions for the cases p=1,...,4 to find all acceptable 
% candidates for a solution z(i)=exp(-b(i)*delta(t)), i=1,...,p, 
% without making any assessments between  bad or good candidates 
% (delta(t)=distance t(j)-t(j-1)). 
% The cases p=1,2 are treated directly by explicit formulas. 
% The cases p=3,4 are treated by a numerical process. 
% The cases p>4 are treated by a heuristic approach. 
%
% Input variables:
% p(1) : length(lambda) asked for. 
%        If length(p)>1 then
%           p(2)=eType, type of exponential sum equation to be used.
%           p(3)=z_max (See output below).  
%           p(4)=sigType, sign combination for equations.
% t, y: a time series (t, y). It must not be equidistant in t.
% S    : Partial sums of y, where y must be eqiudistant in t.
%
% Output variables:
% z    : A vector of solution candidates with entries z(i)=exp(-b(i)*delta(t)).
%        Note that b>0 and real gives 0<z<1. Thus we use a limit z<z_max=1. It
%        can be changed if b>0 is not a demand.
% z_lim: A vector of lower/upper limits for z(i) for p=3 and 4.
%
% Algorithms:
% p=1  : Mean of an ML-estimate of a two-data point estimate.
% p=2-4: A formula based upon equidistant data. Then the exponential
%        sum can be rewritten as a geometrical series. Successive elimination
%        gives an explicit formula for p=2, a one-variable equation for p=3 
%        and a two-variable equation for p=4. 
%        Solutions are found inside the interval [0,1].
% p=3-4: The equations for p=3-4 are a sum of sum(N[j]*z2^j/D[j], j=0,...,4)
%        where N[j] and D[j] are polynomials and 
%        z2 = polynomial + square root of a polynomial. 
%        The denominator D[j] has zeros. 
%        For p=3 they occur as two points which partition the search interval 
%        [0,1] into three parts in which p=3 seem to have exactly one solution. 
%        For p=4 the zeros occur in [0,1]x[0,1] as curves,
%        which follow a certain pattern. 
%        To avoid extra numerical difficulties, the search is made outside 
%        the most complicated part of this pattern. The search is made by 
%        scanning along lines (one fixed variable) so that the situation 
%        becomes like p=3. A difference is that p=4 can have more than 
%        one solution in each subinterval.
% p>4  : For the cases p>4 we use a naïve approach. We construct a
%        p=3-denomintator. Then we partition [0, 1] into three subintervals
%        by the roots of the denominator. We put different numbers of solutions 
%        into different subintervals.
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomopt.com
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written Aug 15, 1997.  Last modified June 21, 1999.
%

function [z, z_lim] = exp_geo(p, t, y)

nargin;
% Or function [z, z_lim] = exp_geo(p, S)

[p,eType,z_max,sigType] = pezs(p);

z_eq=[];
if nargin==3;
   m=length(t);
   eps1=eps;
   if p==1
      % Safe guarded against zeros and changes in sign of y. 
      % These are weighted down to zero.
      index1=1;
      for j=1:m-1
         if sign(y(j+1))*sign(y(j))>0;
            l_est(index1)=y(j+1)/y(j);
            index1=index1+1;
         end
      end
      z=0;
      if index1~=1, z=sum(l_est)/index1; end
   else %p>1
      %Construction of equidistant times.
      ni=length(t);
      delta=t(2)-t(1);
      te=[t(1):delta:t(ni)]';
      ne=length(te);
      %ML-reconstruction of interpolating equidistant data Ye.
      Ye=y(1:2);
      for k=2:ni-1;
         if length(Ye)<ne;
            %This if-construction prevents out of index in te if
            %there is an interval [t(k:k+1)] beyond [te(ne-1:ne)].
            tsubint=[te(length(Ye)+1):delta:t(k+1)]';
            if length(tsubint)>0
               %This if-construction is needed if next te(j) is outside t(k:k+1).
               if y(k)*y(k+1)<0;
                  %This if-construction is needed as the algorithm does not
                  %cope a change of sign in y. In such case use linear interpolation.
                  Ye=[Ye; ( (y(k+1)-y(k))/(t(k+1)-t(k)) ) * ( tsubint-t(k) ) + y(k)];
               else
                  be=log( y(k+1)/y(k) )/( t(k)-t(k+1) );
                  ae=y(k)*exp(be*t(k));
                  Ye=[Ye; ae*exp(-be*tsubint)];
               end
            end
         end   
      end
      
      % Divide the series into partial sums S(i). 
      n=m-rem(m, 2*p);
      q=n/(2*p);
      S=zeros(2*p, 1);
      for i=1:2*p;
         S(i)=sum(y(1+(i-1)*q:i*q));
      end
   end
else %if nargin==2;
   S=t;  %S is stored in t when the call is exp_geo(p, S);
   z = p1z(p, eType, z_max,S);
end
   
if p==2
   z = p2e(eType,S,z_max);
end %if p==2

if p==3 & eType==1
   [Q,P,D,R0,R1,R2,R3,R4] = p3eT1(S);
   % Compute roots D_root(1:2) for the denominator and safe guard 
   % against non real roots and roots outside [0, 1].
   D_root=sort(max(eps, min(real(roots(D)), z_max)));
   %If the roots coincide (due to truncated complex part), separate them a bit. 
   if D_root(1)==D_root(2);
      D_root(2)=D_root(2)+(1-D_root(2))*0.1;
      D_root(1)=D_root(1)*0.9;
   end
   z_lim=D_root;
   sig=[0, 2];
   % If I=[x,x], a divide by zero could occur, but the answer is 
   % a minimum at x as expected.
   sig(1)=any(sigType==[1 2 3 4]);
   I=[0, D_root(1)];
   z(1) = exp_root(I, R0, R1, R2, R3, R4, P, Q, D, sig);
   sig(1)=any(sigType==[1 2 5 6]);
   I=[D_root(1), D_root(2)];
   z(2) = exp_root(I, R0, R1, R2, R3, R4, P, Q, D, sig);
   sig(1)=any(sigType==[1 3 5 7]);
   I=[D_root(2), 1];
   z(3) = exp_root(I, R0, R1, R2, R3, R4, P, Q, D, sig);
end %if p==3 
if p==3 & eType==2

   [Q,P,D,R0,R1,R2,R3,R4] = p3eT2(S);

   % Compute roots D_root(1:2) for the denominator and safe guard 
   % against non real roots and roots outside [0, 1].
   D_root=sort(max(eps, min(real(roots(D)), z_max)));
   %If the roots coincide (due to truncated complex part), separate them a bit. 
   if D_root(1)==D_root(2);
      D_root(2)=D_root(2)+(1-D_root(2))*0.1;
      D_root(1)=D_root(1)*0.9;
   end
   z_lim=D_root;
   sig=[0, 2];
   % If I=[x,x], a divide by zero could occur, but the answer is 
   % a minimum at x as expected.
   I=[eps, D_root(1)];
   z(1) = D_root(1)/2; %First root.
   
   %Find z(2), which is a root with sign change + -> -.
   %Use sig(1)==0 unless testing. It gives better results.
   %
   % If I=[x,x], a divide by zero could occur, but the answer is 
   % a minimum at x as expected.
   
   z_cand = f0min('exp_eq',D_root(1), mean(D_root),...
      [], R0, R1, R2, R3, R4, P, Q, D, [0,2]);         
   zm= exp_eq(0.05*D_root(1)+0.95*z_cand, R0, R1, R2, R3, R4, P, Q, D, sig(1));
   zp= exp_eq(0.05*D_root(2)+0.95*z_cand, R0, R1, R2, R3, R4, P, Q, D, sig(1));
   if (zm>0 & zp<0)
      z(2)=z_cand;
   else
      z_cand = f0min('exp_eq', mean(D_root), D_root(2),...
         [], R0, R1, R2, R3, R4, P, Q, D, [0,2]);         
      zm= exp_eq(0.05*D_root(1)+0.95*z_cand, R0, R1, R2, R3, R4, P, Q, D, sig(1));
      zp= exp_eq(0.05*D_root(2)+0.95*z_cand, R0, R1, R2, R3, R4, P, Q, D, sig(1));
      if (zm>0 & zp<0)
         z(2)=z_cand;
      else
         z_cand = f0min('exp_eq', D_root(1), D_root(2),...
            [], R0, R1, R2, R3, R4, P, Q, D, [0,2]);         
         zm= exp_eq(0.05*D_root(1)+0.95*z_cand, R0, R1, R2, R3, R4, P, Q, D, sig(1));
         zp= exp_eq(0.05*D_root(2)+0.95*z_cand, R0, R1, R2, R3, R4, P, Q, D, sig(1));
         if (zm>0 & zp<0) | (zm*zp<0);
            z(2)=z_cand;
         else            
            z(2)=0.95*D_root(1)+0.05*D_root(1);
         end
      end
   end
   
      
   %if 0 %This algorithm is not as good as the one above.
   %sig(1)=any(sigType==[1 2 5 6]);
   %zint=[D_root(1):(D_root(2)-D_root(1))/50:D_root(2)];      
   %Ysim = exp_eq(zint(2:length(zint)-1), R0, R1, R2, R3, R4, P, Q, D, sig(1));
   %%Do not use the ends z2(1) and z2(length(/z2)) as a division by 0 occur.
   %Ysig=sign(Ysim);
   %%Ysig = sign(exp_eq(z2, R0, R1, R2, R3, R4, P, Q, D, sig(1)));
   %ch_sig=find(Ysig(2:length(Ysig))-Ysig(1:length(Ysig)-1)<0);
   %if length(ch_sig)>0
   %   za=zint(2)+(ch_sig(1)-1)*(zint(length(zint))-zint(1))/50;
   %   zb=zint(2)+ch_sig(1)*(zint(length(zint))-zint(1))/50;
   %   z(2)=f0min('exp_eq',za, zb,...
   %      [], R0, R1, R2, R3, R4, P, Q, D, [1,2]);         
   %else %length(ch_sig)~>0
   %   ch_sig=find(Ysig(2:length(Ysig))-Ysig(1:length(Ysig)-1)>0);
   %   if length(ch_sig)>0
   %      z(2)=0.95*D_root(1)+0.05*D_root(2);
   %   else
   %      z(2)=f0min('exp_eq', D_root(1), D_root(2),...
   %         [], R0, R1, R2, R3, R4, P, Q, D, [sig(1),2]);
   %   end
   %end      
   %end %if 0

   %if 0 %This algorithm is not as good as the first one.
   %I=[D_root(1), D_root(1)+(D_root(2)-D_root(1))/2];
   %z_cand(1,1) = exp_root(I, R0, R1, R2, R3, R4, P, Q, D, sig);
   %z_cand(1,1)= min( max( I(1), z_cand(1, 1)), I(2));
   %z_cand(1,2)= exp_eq(z_cand(1, 1), R0, R1, R2, R3, R4, P, Q, D, sig);
   %I=[D_root(1)+(D_root(2)-D_root(1))/2, D_root(2)];
   %z_cand(2,1) = exp_root(I, R0, R1, R2, R3, R4, P, Q, D, sig);   
   %z_cand(2,1)= min( max( I(1), z_cand(1, 1)), I(2));
   %z_cand(2,2)= exp_eq(z_cand(1, 1), R0, R1, R2, R3, R4, P, Q, D, sig);
   %[trash, j]=min(z_cand(:, 2));    %Find next best root. 
   %z(2)=z_cand(1, 1); %Store next best root in z. 
   %end %if 0

   sig(1)=any(sigType==[1 3 5 7]);
   I=[D_root(2), 1];
   z(3) = exp_root(I, R0, R1, R2, R3, R4, P, Q, D, sig);
   %Third root.
   z=sort(z(:));
end %if p==3 & eType==2

if (p==4) & eType==1 & (length(S)==8);

   [Q,P,D,R0,R1,R2,R3,R4] = p4eT1S8(S);
   
   %We will search zeros z3 and z4 by scanning along lines 
   %0=exp_roots(z3=constant, z4) for some constants stored in the vector z_line.
   %z_line is computed below:
   
   %Eight debug rows: (1998-01-30).
   %Problems occur for problem 13.
   %fprintf('S(2)= %20.18f \n', S(2));
   %fprintf('S(4)= %20.18f \n', S(4));
   %fprintf('S(3)= %25.23f \n', S(3));
   %fprintf('S(2)*S(4)-S(3)^2= %25.23f \n', S(2)*S(4)-S(3)^2);
   %fprintf('D(1,3)= %25.23f \n', D(1,3));
   %Computed by hand it become about -1e-18, truncated to 0. This gives failure.
   %Conclusion; The method is unstable. Try heuristic method instead.
   trash1=roots(D(3, 1:3))';
   trash2=roots(sum(D))';
   if length(trash1)<2 | length(trash2)<2
      %It can happend that this method fails due to D(3,1)==0 which
      %gives only a first order equation with only one root. Then
      %try the heuristic method "if (p>4) | ( (p==4) & (length(S)==6) );"
      %below.
      S=S(1:6);
   else
      D_line=[trash1;...
         zeros(2, 2);...
         trash2];
   %The first and last rows can also be written as roots([0,0,1]*D) and
   %roots([1,1,1]*D);
   %Ignore non real parts and parts outside [0, z_max]:
   D_line=sort(max(eps, min(real(D_line), z_max))')';
   %Now compute the scanning lines z_line and their roots safe guarded
   %from being non real and outside [0, z_max].
   z_line=[0, 0.5*min(D_line(:, 1)), 0.5+0.5*max(D_line(:, 2)), 1];
   D_line(2, :)=roots([z_line(2)^2, z_line(2), 1]*D)';
   D_line(3, :)=roots([z_line(3)^2, z_line(3), 1]*D)';
   D_line=sort(max(eps, min(real(D_line), z_max))')';
   z_lim=[min(D_line(:, 1)), max(D_line(:, 2))];
   z_line=vander(z_line);
   z_line=[z_line(:, 1).*z_line(:, 3), z_line];
   %Now z_line is a matrix with entries z_line(i, j)=z_line(i, 4)^(5-j).
   
   %Now look for solution candidates by finding solutions 
   %to min(log|equation(z3=constant, z4)|) in the 
   %subintervals [D_line(x1, y1), D_line(x1, y2)].
   %Not every interval is searched through but only:
   % Two best roots in [0, D_line(j, 1)] where j=1:4,
   % Two best roots in [D_line(j, 1), D_line(j, 2)] where j=2:4,
   % and among these four roots (later) pick out the three best.
   % One root in [D_line(4, 2), z_max].
   %The results are stored in a matrix z_cand1 and best roots sorted into z_cand2.
   sig=[0, 2];  %Initiate a parameter for exp_root.
   
   %Find two best roots of four in [0, D_line(j, 1)].
   z_cand1=[];

   sig(1)=any(sigType==[1 2 3 4]);   
   for i=1:4
      % Determine the other coefficient matrices.
      P4=P*(z_line(i, 3:5))';
      D4=P*(z_line(i, 3:5))';
      Q4=Q*(z_line(i, :))';
      R04=R0*(z_line(i, :))';
      R14=R1*(z_line(i, :))';
      R24=R2*(z_line(i, :))';
      R34=R3*(z_line(i, 3:5))';
      R44=R4*(z_line(i, 4:5))';
      % If I=[x,x], a divide by zero could occur, but the answer is 
      % a minimum at x as expected.
      I=[eps, D_line(i, 1)];
      z_cand1(i, 1)= exp_root(I, R04, R14, R24, R34, R44, P4, Q4, D4, sig);
      z_cand1(i, 1)= min( max( I(1), z_cand1(i, 1)), I(2));
      z_cand1(i, 2)= exp_eq(z_cand1(i, 1), R04, R14, R24, R34, R44,...
         P4, Q4, D4, sig);
   end
   [trash, j]=min(z_cand1(:, 2));    %Find best root. 
   z_cand2(1,1:2)=z_cand1(j, :); %Store best root in z_cand2. 
   z_cand1(j,2)=inf;             %Destroy best root in z_cand1.  
   [trash, j]=min(z_cand1(:, 2));    %Find next best root. 
   z_cand2(2,1:2)=z_cand1(j, :); %Store next best root in z_cand2. 
   %New lagorithm: take one root in each subinterval and one extra fast 
   %and well separated root. Store it in z_sol.
   if min(z_cand2(:,1))>0.5*max(z_cand2(:,1));
      z_sol=[0.5*z_cand2(1,1);z_cand2(1,1)];
   else
      z_sol=[z_cand2(1:2,1)];
   end
      
   %Find two best roots of three in [D_line(j, 1), D_line(j, 2)].
   z_cand1=[];
   sig(1)=any(sigType==[1 2 5 6]);   
   for i=2:4
      % Determine the other coefficient matrices.
      P4=P*(z_line(i, 3:5))';
      D4=P*(z_line(i, 3:5))';
      Q4=Q*(z_line(i, :))';
      R04=R0*(z_line(i, :))';
      R14=R1*(z_line(i, :))';
      R24=R2*(z_line(i, :))';
      R34=R3*(z_line(i, 3:5))';
      R44=R4*(z_line(i, 4:5))';
      % If I=[x,x], a divide by zero could occur, but the answer is 
      % a minimum at x as expected.
      I=[D_line(i, 1), D_line(i, 2)];
      z_cand1(i-1, 1) = exp_root(I, R04, R14, R24, R34, R44, P4, Q4, D4, sig);
      z_cand1(i-1, 1)= min( max( I(1), z_cand1(i-1, 1)), I(2));
      z_cand1(i-1, 2) = exp_eq(z_cand1(i-1, 1), R04, R14, R24, R34, R44,...
         P4, Q4, D4, sig);
      
   end
   [trash, j]=min(z_cand1(:, 2));    %Find best root. 
   z_cand2(3,1:2)=z_cand1(j, :); %Store best root in z_cand2. 
   z_cand1(j,2)=inf;             %Destroy best root in z_cand1.  
   [trash, j]=min(z_cand1(:, 2));    %Find next best root. 
   z_cand2(4,1:2)=z_cand1(j, :); %Store next best root in z_cand2. 
   
   z_sol=[z_sol; z_cand2(3,1)];
         
   %Pick out three best roots of four.
   [trash, j]=min(z_cand2(:, 2)); %Find best root and store it in z. 
   z(1)=z_cand2(j, 1);
   z_cand2(j,2)=inf;             %Destroy best root in z_cand1.  
   [trash, j]=min(z_cand2(:, 2)); %Find next best root and store it in z. 
   z(2)=z_cand2(j, 1);
   z_cand2(j,2)=inf;             %Destroy next best root in z_cand1.  
   [trash, j]=min(z_cand2(:, 2)); %Find next next best root and store it in z. 
   z(3)=z_cand2(j, 1);
   
   
   %Find root in [D_line(j, 2), z_max].
   z_cand1=[];
   sig(1)=any(sigType==[1 3 5 7]);
   i=4;
   % Determine the other coefficient matrices.
   P4=P*(z_line(i, 3:5))';
   D4=P*(z_line(i, 3:5))';
   Q4=Q*(z_line(i, :))';
   R04=R0*(z_line(i, :))';
   R14=R1*(z_line(i, :))';
   R24=R2*(z_line(i, :))';
   R34=R3*(z_line(i, 3:5))';
   R44=R4*(z_line(i, 4:5))';
   % If I=[x,x], a divide by zero could occur, but the answer is 
   % a minimum at x as expected.
   I=[D_line(i, 2), z_max];
   z(i)= exp_root(I, R04, R14, R24, R34, R44, P4, Q4, D4, sig);
   z(i)= min( max( I(1), z(i)), I(2));
   z_sol=[z_sol; z(i)];

   % Now we have four solution candidates
   % in two versions. One is z and one is z_sol. 
   %z_sol is the better one as it guarantees three good roots
   %and a well separated fourth root.
   z=z_sol(:);
   end %if length(trash1)<2 | length(trash2)<2 

end % if (p==4)...

if (p>4) | ( (p==4) & (length(S)==6)) | (eType==2 & p>2);
   Pmat=[floor(0.4*p), 1];
   Pmat=[p-sum(Pmat), Pmat];
   if eType==1;
      D = [S(2)^2-S(3)*S(1), S(4)*S(1)-S(3)*S(2), S(3)^2-S(4)*S(2)];
   elseif eType==2;
      D = [S(2)^2-S(2)*S(3)-S(2)*S(4)+S(4)*S(1)+S(3)^2-S(1)*S(3),...
      -S(2)*S(3)+S(2)*S(5)-S(5)*S(1)-S(3)*S(4)+S(4)*S(1)+S(3)^2,...
      S(2)*S(5)-S(2)*S(4)+S(3)^2-S(3)*S(4)-S(5)*S(3)+S(4)^2];
   end
   % Compute roots D_root(1:2) for the denominator and safe guard 
   % against non real roots and roots outside [0, 1].
   D_root=sort(max(eps, min(real(roots(D)), z_max)));
   %If the roots coincide (due to truncated complex part), separate them a bit. 
   if D_root(1)==D_root(2);
      D_root(2)=D_root(2)+(1-D_root(2))*0.1;
      D_root(1)=D_root(1)*0.9;
   end
   %Compute heuristic estimate of z.
   z(1:Pmat(1))=D_root(1)/(Pmat(1)+1)*[1:Pmat(1)];
   z(Pmat(1)+1:p-1)=D_root(1)+(D_root(2)-D_root(1))/(Pmat(2)+1)*[1:Pmat(2)];
   z(p)=D_root(2)+0.12*(z_max-D_root(2));
end %if (p>4)...
if ~exist('z','var');
   z=[1/(p+1):1/(p+1):1-1/(p+1)];
end
z=sort(z(:));
if ~exist('z_lim','var');
   z_lim=[];
end
z_lim=z_lim(:);

% end % of routine exp_geo
% MODIFICATION LOG
%
% 981022  hkh  Removed global n_f statement

% =======================================================================
%        exp_root.m
% ========================================================================
%
% function z = exp_root(I, R0, R1, R2, R3, R4, P, Q, D, sig);
%
% The function tries to find a root of an equation of the form 
% 
%     0 = N^4*R4/(2D^3) + N^3*R3/(2D^2) + N^2*R2/D^2 + N*R1/D) + R0, 
% where N=P+sqrt(Q) and R0,...,D are polynomials in z of degree less than four.
%
% Input variables
% I        : An interval [a,b] to search. 
%            (Between asymptotes from the denominator D).
% R0,...,D : Matrices with polynomial coefficients
% sig      : A sign which controls whether to use N = P + or - sqrt(Q).
% z        : A root.
%
% Algorithms
% As the equation is fractional with polynonial coefficients and a square root,
% the derivatives are complicated. 
% Thus a derivative free method seems to be the best choice.
% The application for this equation (finding initial values for the 
% fitting of exponential sums to empirical data) has very small numerical 
% values and also work in very different scales. 
% About 10^{-5} downto 10^{-20} for different problems is normal
% with an interval of five to ten powers for each certain root search. 
% Taking the logarithm y=log(|equation|) of the equation will transform the 
% problem to a minimization problem. 
% In many cases the curve bends slowly over the interval and has a sharper 
% "arrow" downwards close to a zero.
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomopt.com
% Copyright (c) 1997-2004 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written May 26, 1997.      Last modified Dec 5, 2001.
%

%Further developement.
% In the case p=3 (see exp_geo.m) there seems to be only one root inside each
% interval. For p=4 there could be more than one root. Added a special procedure
% for handling this!
%

function z = exp_root(I, R0, R1, R2, R3, R4, P, Q, D, sig, eType)

if nargin<11; eType=1; end; %eType=2 is treated the same way.

if eType<3
   z = f0min('exp_eq',I(1), I(2), [], R0, R1, R2, R3, R4, P, Q, D, sig);
elseif eType>2 & eType<5
   P_root=I(1); z_max=I(2);
   %Find z(1), which is a root with sign change - -> +.
   Prob.R0  = R0;
   Prob.R1  = R1;
   Prob.R2  = R2;
   Prob.R3  = R3;
   Prob.R4  = R4;
   Prob.P   = P;
   Prob.Q   = Q;
   Prob.D   = D;
   Prob.sig = sig;
   Prob.FUNCS.f0 = 'exp_eq0';

   %fz = exp_eq(0, R0, R1, R2, R3, R4, P, Q, D, [1]);
   fz = exp_eq0(0, Prob);

   z  = Tfzero(-20,20,Prob,0);

   %if findstr('5.3',version)
   %   z=fzero('exp_eq',0, [], R0, R1, R2, R3, R4, P, Q, D, [1]);
   %else
   %   z=fzero('exp_eq',0, [], [], R0, R1, R2, R3, R4, P, Q, D, [1]);
   %end

      if fz<0;
         if z<0 | P_root<z; %Fail in search.
            z=f0min('exp_eq',0, P_root*0.95,...
               [], R0, R1, R2, R3, R4, P, Q, D, [1,2]);
         end
      else %fz>0.
         %if 0<z & z<P_root; %If success in search.
         %   I2=[z+(P_root-z)*0.05, z+(P_root-z)*0.95];
         %else %If failure in search, e.g. z=NaN.
            z=f0min('exp_eq',0, P_root*0.95, [],...
               R0, R1, R2, R3, R4, P, Q, D, [1,2]);
            I2=[z+(P_root-z)*0.05, z+(P_root-z)*0.95];
         %end
         z=f0min('exp_eq',I2(1), I2(2), [], R0, R1, R2, R3, R4, P, Q, D, [1,2]);
      end
      za=z(1)*0.98;
      zb=z(1)*0.98+z_max*0.02;
      fza = exp_eq(za, R0, R1, R2, R3, R4, P, Q, D, [1]);
      fzb = exp_eq(zb, R0, R1, R2, R3, R4, P, Q, D, [1]);
      if ~(fza<0 & fzb>0)
         z=f0min('exp_eq',zb, P_root, [], R0, R1, R2, R3, R4, P, Q, D, [1,2]);
      end
      %z(1) found.
      
      %Find z(2), which is a root with sign change + -> -.
      z1=0.95*z(1)+0.05*z_max;      
      z2=[z1:(z_max-z1)/50:z_max];
      Ysim = exp_eq(z2, R0, R1, R2, R3, R4, P, Q, D, [1]);
      Ysig=sign(Ysim);
      %Ysig = sign(exp_eq(z2, R0, R1, R2, R3, R4, P, Q, D, [1]));
      ch_sig=find(Ysig(2:length(Ysig))-Ysig(1:length(Ysig)-1)<0);
      if length(ch_sig)==1
         za=z1+(ch_sig-1)*(z_max-z1)/50;
         zb=z1+ch_sig*(z_max-z1)/50;
         z(2)=f0min('exp_eq',za, zb,...
            [], R0, R1, R2, R3, R4, P, Q, D, [1,2]);         
      else %length(ch_sig)~=1
         z(2)=f0min('exp_eq',z1, z_max,...
            [], R0, R1, R2, R3, R4, P, Q, D, [1,2]);                        
      end      
      
   
   
   %if 0
   %P_root=I(1);
   %z_max=I(2);
   %%Find z(1), which is a root with sign change - -> +.
   %fz = exp_eq(0, R0, R1, R2, R3, R4, P, Q, D, [1]);
   %z=fzero('exp_eq',I(1), [], [], R0, R1, R2, R3, R4, P, Q, D, [1]);
   %if fz<0;
   %   if z<0 | P_root<z; %Fail in search.
   %      z=f0min('exp_eq',0, P_root*0.95,...
   %         [], R0, R1, R2, R3, R4, P, Q, D, [1,2]);
   %   end
   %else %fz>0.
   %   if 0<z & z<P_root; %If success in search.
   %      I2=[z+(P_root-z)*0.05, z+(P_root-z)*0.95];
   %   else %If failure in search, e.g. z=NaN.
   %      z=f0min('exp_eq',I(1), P_root*0.95, [],...
   %         R0, R1, R2, R3, R4, P, Q, D, [1,2]);
   %      I2=[z+(P_root-z)*0.05, z+(P_root-z)*0.95];
   %   end
   %   z=f0min('exp_eq',I2(1), I2(2), [], R0, R1, R2, R3, R4, P, Q, D, [1,2]);
   %end
   %za=z(1)*0.98;
   %zb=z(1)*0.98+z_max*0.02;
   %fza = exp_eq(za, R0, R1, R2, R3, R4, P, Q, D, [1]);
   %fzb = exp_eq(zb, R0, R1, R2, R3, R4, P, Q, D, [1]);
   %if ~(fza<0 & fzb>0)
   %   z=f0min('exp_eq',zb, P_root, [], R0, R1, R2, R3, R4, P, Q, D, [1,2]);
   %end
   %%z(1) found.
  % 
  % %Find z(2), which is a root with sign change + -> -.
  % z1=0.95*z(1)+0.05*z_max;      
  % z2=[z1:(z_max-z1)/50:z_max];
  % Ysim = exp_eq(z2, R0, R1, R2, R3, R4, P, Q, D, [1]);
  % Ysig=sign(Ysim);
  % %Ysig = sign(exp_eq(z2, R0, R1, R2, R3, R4, P, Q, D, [1]));
  % ch_sig=find(Ysig(2:length(Ysig))-Ysig(1:length(Ysig)-1)<0);
  % if length(ch_sig)==1
  %    za=z1+(ch_sig-1)*(z_max-z1)/50;
  %    zb=z1+ch_sig*(z_max-z1)/50;
  %    z(2)=f0min('exp_eq',za, zb,...
  %       [], R0, R1, R2, R3, R4, P, Q, D, [1,2]);         
  % else %length(ch_sig)~=1
  %    z(2)=f0min('exp_eq',z1, z_max,...
  %       [], R0, R1, R2, R3, R4, P, Q, D, [1,2]);                        
  % end 
  % end %if 0
end 

% =================================================
function [Q,P,D,R0,R1,R2,R3,R4] = p3eT2(S)
% =================================================

   Q=Qcomp(S);

   P = [S(3)*S(4)-S(2)*S(5)+S(5)*S(1)+S(2)*S(3)-S(3)^2-S(4)*S(1),...
         2*S(3)*S(4)-S(4)^2-S(2)*S(5)+S(2)*S(6)-S(6)*S(1)+S(5)*S(1)-...
         S(3)^2,...
         S(5)*S(4)-S(2)*S(5)+S(2)*S(6)-S(6)*S(3)+S(3)*S(4)-S(4)^2];
   D = 2*[S(2)^2-S(2)*S(3)-S(2)*S(4)+S(4)*S(1)+S(3)^2-S(1)*S(3),...
         -S(2)*S(3)+S(2)*S(5)-S(5)*S(1)-S(3)*S(4)+S(4)*S(1)+S(3)^2,...
         S(2)*S(5)-S(2)*S(4)+S(3)^2-S(3)*S(4)-S(5)*S(3)+S(4)^2];
   
   R0 =[-S(4)*S(2)^2+S(2)^2*S(5)-S(3)^3+S(3)^2*S(2)+S(4)*S(3)^2-2*...
         S(3)*S(5)*S(2)+S(4)^2*S(2)+S(5)*S(3)^2-S(3)*S(4)^2,...
         -S(4)*S(2)^2+S(2)^2*S(5)-S(3)^3+S(3)^2*S(2)+S(4)*S(3)^2-2*...
         S(3)*S(5)*S(2)+S(4)^2*S(2)+S(5)*S(3)^2-S(3)*S(4)^2,...
         -2*S(3)*S(5)*S(2)-2*S(2)*S(4)*S(5)+3*S(5)*S(3)^2-2*S(3)*S(5)^...
         S(4)-S(2)^2*S(7)+S(2)^2*S(5)+2*S(5)*S(4)^2-S(3)^3+S(4)*...
         S(3)^2-S(4)^3+S(4)^2*S(2)+2*S(2)*S(7)*S(3)+S(2)*S(5)^2-S(7)*...
         S(3)^2-S(3)*S(5)^2,...
         S(2)*S(5)^2+S(4)^2*S(2)-2*S(3)*S(5)*S(4)+S(3)*S(5)^2+3*S(5)*...
         S(4)^2-2*S(3)*S(5)*S(2)+S(5)*S(3)^2+2*S(2)*S(7)*S(3)-2*S(4)*...
         S(5)^2+S(4)*S(3)^2-S(3)*S(4)^2-S(4)^3-2*S(7)*S(3)^2-2*S(2)*...
         S(7)*S(4)+2*S(7)*S(4)*S(3),...
         -S(7)*S(3)^2-S(3)*S(4)^2-2*S(4)*S(5)^2+2*S(5)*S(4)^2+S(5)*...
         S(3)^2+S(5)^3+2*S(7)*S(4)*S(3)-S(3)*S(5)^2-S(7)*S(4)^2];

R1=R1comp(S);

R2(1)= -3*S(3)*S(2)^2-S(4)*S(2)^2+S(2)^2*S(5)-S(3)^3+3*S(3)^2*...
      S(2)+S(1)^2*S(5)-S(4)*S(1)^2+2*S(4)*S(2)*S(1)-2*S(2)*S(5)*S(1)+...
      S(2)^3;

R2(2)=-S(1)*S(3)^2-S(4)^2*S(1)-2*S(3)*S(5)*S(2)+3*S(2)^2*S(5)-2*...
      S(3)^3+6*S(3)^2*S(2)+S(4)*S(3)^2-2*S(4)*S(2)^2+S(4)^2*S(2)-4*...
      S(3)*S(2)^2-S(4)*S(1)^2+S(1)^2*S(5)-4*S(2)*S(5)*S(1)+2*S(3)*...
      S(5)*S(1)-2*S(4)*S(3)*S(2)+S(2)^3+4*S(4)*S(2)*S(1);

R2(3)=6*S(3)^2*S(2)-4*S(2)*S(5)*S(1)+S(1)^2*S(5)-6*S(3)*S(5)*...
      S(2)-2*S(2)*S(4)*S(5)+3*S(5)*S(3)^2-S(2)^2*S(7)+6*S(2)^2*...
      S(5)-2*S(4)^2*S(1)-3*S(3)*S(2)^2-3*S(3)^3+2*S(4)^2*S(2)-2*S(4)*...
      S(2)^2-S(7)*S(1)^2+2*S(5)*S(4)*S(1)+2*S(2)*S(7)*S(1)+...
      2*S(4)*S(2)*S(1);

R2(4)=3*S(3)^2*S(2)+3*S(2)^2*S(5)-2*S(2)*S(5)*S(1)+2*S(5)*S(4)*S(1)+...
      S(2)*S(5)^2+2*S(4)^2*S(2)-2*S(3)*S(5)*S(4)-6*S(3)*S(5)*S(2)+3*...
      S(5)*S(3)^2+2*S(2)*S(7)*S(3)+2*S(3)*S(4)^2-S(4)^3-S(4)*S(2)^2-...
      S(1)*S(5)^2-2*S(2)^2*S(7)-S(4)^2*S(1)-2*S(4)*S(3)*S(2)+2*S(3)*...
      S(5)*S(1)+2*S(2)*S(7)*S(1)-2*S(7)*S(3)*S(1)-2*S(3)^3;

R2(5)=S(2)^2*S(5)-S(2)^2*S(7)-2*S(2)*S(4)*S(5)-S(7)*S(3)^2+2*S(2)*...
      S(7)*S(3)+S(2)*S(5)^2-2*S(3)*S(5)*S(4)+2*S(5)*S(4)^2+3*S(5)*...
      S(3)^2-S(3)^3-2*S(3)*S(5)*S(2)+S(4)^2*S(2)-S(4)^3-S(3)*S(5)^2+...
      S(4)*S(3)^2;

   R3 = [-S(3)+S(2), 2*S(2)-S(1)-S(3), S(2)-S(1)]/2;
   R4 = [ -S(3)+S(2),  S(2)-S(1)]/2;

% =================================================
function [Q,P,D,R0,R1,R2,R3,R4] = p4eT1S8(S)
% =================================================
   % Determine coefficients in the equation. 
   % The coefficient matrices are symmetric.
   % Compute only the upper half and then set the lower part by symmetry.
   P = [S(1)*S(4)-S(2)*S(3), S(3)^2-S(1)*S(5), S(2)*S(5)-S(4)*S(3);...
         0, S(2)*S(5)+S(6)*S(1)-2*S(4)*S(3), S(4)^2-S(6)*S(2);...
         0, 0, S(6)*S(3)-S(5)*S(4)];
   sz=size(P);
   for i=2:sz(1);
      for j=1:i-1;
         P(i,j)=P(j,i);
      end
   end
   D = 2*[S(1)*S(3)-S(2)^2, S(2)*S(3)-S(1)*S(4), double(S(2))*...
         double(S(4))-double(S(3))^2;...
         0, -S(3)^2+S(1)*S(5), -S(2)*S(5)+S(4)*S(3);...
         0, 0, -S(4)^2+S(3)*S(5)];
   sz=size(D);  
   for i=2:sz(1);
      for j=1:i-1;
         D(i,j)=D(j,i);
      end
   end

   Q=Qcomp2(S);


   sz=size(Q);

   for i=2:sz(1);
      for j=1:i-1;
         Q(i,j)=Q(j,i);
      end
   end
   R0= [0,...
        -S(2)*(S(2)*S(4)-S(3)^2),...
         S(2)^2*S(5)-S(3)^3,...
         -2*S(3)*S(5)*S(2)+S(3)^2*S(4)+S(2)*S(4)^2,...
         S(3)*(-S(4)^2+S(3)*S(5));...
       0,...
         S(2)^2*S(5)-S(3)^3,...
         S(3)^2*S(4)-2*S(2)*S(5)*S(3)+S(2)*S(4)^2,...
         S(3)*(S(3)*S(5)-S(4)^2),...
         0;...
       0, 0,...
         -S(7)*S(2)^2+3*S(3)^2*S(5)-2*S(2)*S(5)*S(4),...
         2*S(7)*S(2)*S(3)-S(4)^3-2*S(5)*S(4)*S(3)+S(2)*S(5)^2,...
         -S(7)*S(3)^2+2*S(5)*S(4)^2-S(3)*S(5)^2;...
       0, 0, 0,...
         -2*S(7)*S(4)*S(2)-2*S(7)*S(3)^2+S(3)*S(5)^2+3*S(5)*S(4)^2,...
         2*S(4)*(S(7)*S(3)-S(5)^2);...
       0, 0, 0, 0,...
         -S(7)*S(4)^2+S(5)^3];
   sz=size(R0);
   for i=2:sz(1);
      for j=1:i-1;
         R0(i,j)=R0(j,i);
      end
   end
   R1= [0,...
         2*S(2)*S(4)*S(1)-S(1)*S(3)^2-S(2)^2*S(3),...
         -S(2)*(S(2)*S(4)-3*S(3)^2+2*S(1)*S(5)),...
         2*S(3)*S(5)*S(1)+2*S(2)^2*S(5)-S(3)^3-2*S(3)*S(2)*S(4)-S(1)*S(4)^2,...
         S(2)*S(4)^2+S(3)^2*S(4)-2*S(2)*S(5)*S(3);...
       0,...
         -2*S(2)*(S(2)*S(4)+S(1)*S(5)-2*S(3)^2),...
         -S(1)*S(4)^2+3*S(2)^2*S(5)-2*S(3)*S(2)*S(4)-2*S(3)^3+2*S(3)*S(5)*S(1),...
         -4*S(2)*S(5)*S(3)+2*S(2)*S(4)^2+2*S(3)^2*S(4),...
         S(3)*(S(3)*S(5)-S(4)^2);...
       0, 0,...
         2*S(2)*S(4)^2+2*S(7)*S(2)*S(1)-6*S(2)*S(5)*S(3)+2*S(1)*S(5)*S(4),...
         -2*S(7)*S(2)^2-S(1)*S(5)^2+3*S(3)^2*S(5)+2*S(3)*S(4)^2-2*S(7)*S(3)*S(1),...
         -S(4)^3+S(2)*S(5)^2-2*S(5)*S(3)*S(4)+2*S(7)*S(3)*S(2);...
       0, 0, 0,...
         -6*S(5)*S(3)*S(4)-2*S(4)^3+6*S(7)*S(3)*S(2)+2*S(7)*S(4)*S(1),...
         -2*S(7)*S(3)^2+3*S(5)*S(4)^2+S(3)*S(5)^2-2*S(7)*S(4)*S(2);...
       0, 0, 0, 0,...
         2*S(4)*(-S(5)^2+S(7)*S(3))];
   sz=size(R1);
   for i=2:sz(1);
      for j=1:i-1;
         R1(i,j)=R1(j,i);
      end
   end
   R2= [0,...
         S(2)^3-S(1)^2*S(4),...
         S(1)^2*S(5)-3*S(3)*S(2)^2+2*S(1)*S(2)*S(4),...
         -S(2)*(S(2)*S(4)+2*S(1)*S(5)-3*S(3)^2),...
         S(2)^2*S(5)-S(3)^3;...
       0,...
         S(1)^2*S(5)-4*S(3)*S(2)^2-S(1)*S(3)^2+4*S(1)*S(2)*S(4),...
         -2*S(2)*(S(2)*S(4)+2*S(1)*S(5)-3*S(3)^2),...
         3*S(2)^2*S(5)-2*S(3)^3-2*S(3)*S(2)*S(4)-S(1)*S(4)^2+2*S(1)*S(5)*S(3),...
         -2*S(2)*S(5)*S(3)+S(2)*S(4)^2+S(3)^2*S(4);...
       0, 0,...
         -S(7)*S(1)^2-3*S(3)^3+6*S(2)^2*S(5)-2*S(1)*S(4)^2,...
         2*S(7)*S(1)*S(2)+2*S(2)*S(4)^2-6*S(2)*S(5)*S(3)+2*S(1)*S(5)*S(4),...
         -2*S(2)*S(5)*S(4)-S(7)*S(2)^2+3*S(5)*S(3)^2;...
       0, 0, 0,...
         -S(1)*S(5)^2-2*S(7)*S(1)*S(3)-2*S(7)*S(2)^2+2*S(3)*S(4)^2+3*S(5)*S(3)^2,...
         S(2)*S(5)^2+2*S(7)*S(3)*S(2)-S(4)^3-2*S(5)*S(3)*S(4);...
       0, 0, 0, 0,...
         -S(7)*S(3)^2+2*S(5)*S(4)^2-S(3)*S(5)^2];
   sz=size(R2);
   for i=2:sz(1);
      for j=1:i-1;
         R2(i,j)=R2(j,i);
      end
   end
   R3= [0, S(1), -S(2); S(1), -2*S(2), S(3); -S(2), S(3), 0];
   R4= [S(1), -S(2); -S(2), S(3)];

% =================================================
function [Q,P,D,R0,R1,R2,R3,R4] = p3eT1(S)
% =================================================
   % Determine coefficients in the equation. 
   P = [S(2)*S(3)-S(4)*S(1), S(1)*S(5)-S(3)^2, S(4)*S(3)-S(5)*S(2)];
   D = 2*[S(2)^2-S(3)*S(1), S(4)*S(1)-S(3)*S(2), S(3)^2-S(4)*S(2)];
   Q = [S(1)^2*S(4)^2+4*S(2)^3*S(4)-3*S(3)^2*S(2)^2+4*S(3)^3*S(1)-...
         6*S(1)*S(2)*S(3)*S(4),...
         6*S(1)*S(5)*S(2)*S(3)-4*S(2)^3*S(5)+4*S(1)*S(4)^2*S(2)-...
         2*S(1)^2*S(4)*S(5)+2*S(3)^3*S(2)-6*S(3)^2*S(4)*S(1),...
         6*S(3)^2*S(4)*S(2)-3*S(3)^4+6*S(2)^2*S(5)*S(3)-8*S(2)^2*S(4)^2+...
         6*S(3)*S(4)^2*S(1)+S(1)^2*S(5)^2-2*S(2)*S(4)*S(5)*S(1)-...
         6*S(1)*S(5)*S(3)^2,...
         -6*S(2)*S(5)*S(3)^2+6*S(3)*S(4)*S(5)*S(1)+4*S(2)^2*S(4)*S(5)+...
         2*S(3)^3*S(4)-2*S(1)*S(5)^2*S(2)-4*S(1)*S(4)^3,...
         4*S(3)^3*S(5)+4*S(2)*S(4)^3+S(2)^2*S(5)^2-6*S(3)*S(4)*S(5)*S(2)-3*...
         S(3)^2*S(4)^2];
   R0= [S(4)*S(2)^2-S(3)^2*S(2), 0, 2*S(3)^2*S(4)-S(2)*S(4)^2-S(6)*S(2)^2,...
         -2*S(3)*S(4)^2+2*S(6)*S(3)*S(2), -S(6)*S(3)^2+S(4)^3];
   R1= [2*S(1)*S(4)*S(2)-S(2)^2*S(3)-S(1)*S(3)^2, -S(4)*S(2)^2+S(3)^2*S(2),...
         -S(1)*S(4)^2-2*S(6)*S(1)*S(2)+S(3)^3+2*S(2)*S(4)*S(3),...
         -3*S(3)^2*S(4)+2*S(6)*S(3)*S(1)+2*S(6)*S(2)^2-S(2)*S(4)^2,...
         2*S(3)*S(4)^2-2*S(6)*S(3)*S(2)];
   R2= [-S(2)^3+S(1)^2*S(4), -2*S(1)*S(4)*S(2)+S(2)^2*S(3)+S(1)*S(3)^2,...
         -S(6)*S(1)^2+3*S(4)*S(2)^2-2*S(4)*S(1)*S(3),...
         -2*S(2)*S(4)*S(3)+S(1)*S(4)^2+2*S(6)*S(1)*S(2)-S(3)^3,...
         2*S(3)^2*S(4)-S(2)*S(4)^2-S(6)*S(2)^2];
   R3= [S(1), -S(2), 0]/2; 
   R4= [-S(1), S(2)]/2;

% =================================================
function [Q,P,D,R0,R1,R2,R3,R4] = eT3eT4(S)
% =================================================
      P=[S(2), -2*S(3), S(4)];
      D=[S(1), -2*S(2), S(3)];
      Q=[S(2)^2-S(1)*S(3), -2*S(2)*S(3)+2*S(1)*S(4),...
            -2*S(4)*S(2)+3*S(3)^2-S(1)*S(5), -2*S(4)*S(3)+2*S(2)*S(5),...
            S(4)^2-S(3)*S(5)];
      R0=[2*S(3), -3*S(4), 0, S(6)];
      R1=-4*[S(2), -2*S(3), S(4), 0];
      R2=[2*S(1), -7*S(2), 8*S(3), -3*S(4)];
      R3=[2];
      R4=[0];
% =================================================
function z = p1z(p, eType, z_max,S)
% =================================================

   z=[];
   if p==1
      if eType==1
         z=max(eps, min(S(2)/S(1), z_max));
      elseif eType==2
         z=max(eps, min( (S(3)-S(2))/(S(2)-S(1)), z_max));
      elseif eType==3 | eType==4
         if S(1)==0
            z=S(3)/2/S(2);
         else
            z=max(eps, min( real([S(2)-sqrt(S(2)^2-S(1)*S(3))]/S(1)), z_max));
         end
      end         
   end
% =================================================
function [p,eType,z_max,sigType] = pezs(v)
% =================================================

p=1; eType=1; z_max=1; sigType=1;

if length(v) > 0, p=v(1); end
if length(v) > 1, eType=v(2); end
if length(v) > 2, z_max=v(3); end
if length(v) > 3, sigType=v(4); end

% =================================================
function z = p2e(eType,S,z_max)
% =================================================

   if eType==1
      %Solve equations by Maple-formulas based upon geometrical sum formulation.
      %It is a second order eq. with solutions (P +\- sqrt(Q))/D.
      P = double(S(2)*S(3))-double(S(1)*S(4));
      D = 2*( double(S(2))^2-double(S(1))*double(S(3)));
      Q = S(1)^2*S(4)^2+4*S(2)^3*S(4)-3*S(3)^2*S(2)^2+4*S(3)^3*S(1)-...
         6*S(1)*S(2)*S(3)*S(4);
      % Compute solutions.
      sqrtQD=sqrt(max(0, Q))/D;      PD=P/D; 
      z=[PD-sqrtQD;PD+sqrtQD];
      %Safeguard so that the solution is inside [0, 1].
      z=sort(max(eps, min(z, z_max)));   
      %If the roots coincide (due to removed complex parts), separate them a bit.
      if z(1)==z(2);
         z(2)=z(2)+(1-z(2))*0.1;
         z(1)=z(1)*0.9;
      end
   elseif eType==2
      %Solve equations by Maple-formulas based upon geometrical sum formulation.
      %It is a second order eq. with solutions (P +\- sqrt(Q))/D.
      Q = -2*S(1)*S(5)^2*S(2)+4*S(1)*S(4)^2*S(2)+6*S(3)*S(2)^2*S(5)...
         -6*S(3)^2*S(5)*S(2)-6*S(1)*S(5)*S(3)^2+4*S(5)*S(2)^2*S(4)+6*S(1)*...
         S(4)^2*S(3)+6*S(3)^2*S(4)*S(2)-2*S(1)^2*S(4)*S(5)-6*S(1)*S(4)*...
         S(3)^2-3*S(3)^2*S(4)^2+2*S(3)^3*S(4)+4*S(3)^3*S(5)+2*S(3)^3*...
         S(2)-4*S(1)*S(4)^3-3*S(3)^2*S(2)^2+S(1)^2*S(4)^2+S(1)^2*...
         S(5)^2+4*S(3)^3*S(1)-8*S(2)^2*S(4)^2-4*S(2)^3*S(5)+4*S(2)^3*...
         S(4)+4*S(4)^3*S(2)+6*S(1)*S(4)*S(5)*S(3)-2*S(1)*S(4)*S(5)*S(2)+...
         S(5)^2*S(2)^2-6*S(3)*S(4)*S(5)*S(2)-3*S(3)^4-6*S(3)*S(2)*S(1)*...
         S(4)+6*S(3)*S(2)*S(1)*S(5);
      D = 2*(S(3)*S(1)-S(1)*S(4)-S(2)^2+S(3)*S(2)+S(4)*S(2)-S(3)^2);
      P = -S(3)*S(4)+S(3)^2-S(3)*S(2)+S(1)*S(4)+S(5)*S(2)-S(1)*S(5);
      % Compute solutions.
      sqrtQD=sqrt(max(0, Q))/D;      PD=P/D;  
      z=[PD-sqrtQD;PD+sqrtQD];
      %Safeguard so that the solution is inside [0, 1].
      z=sort(max(eps, min(z, z_max)));   
      %If the roots coincide (due to removed complex parts), separate them a bit.
      if z(1)==z(2);
         z(2)=z(2)+(1-z(2))*0.1;
         z(1)=z(1)*0.9;
      end
   elseif eType==3 | eType==4
     [Q,P,D,R0,R1,R2,R3,R4] = eT3eT4(S);
      % Compute roots D_root and P_root and safe guard 
      % against non real roots and roots outside [0, 1].
      D_root=min(real(roots(D)));    P_root=min(real(roots(P)));   
      D_root=max(eps, min(D_root, z_max));
      P_root=max(eps, min(P_root, z_max));
      I=[P_root, z_max];
      sig=[1];
      % If I=[x,x], a divide by zero could occur, but the answer is 
      % a minimum at x as expected.
      z = exp_root(I, R0, R1, R2, R3, R4, P, Q, D, sig, eType);
   end %if eType...

% =================================================
function Q=Qcomp(S)
% =================================================
Q(1)=Qcomp11(S);
Q(2)=Qcomp12(S);
Q(3)=Qcomp13(S);
Q(4)=Qcomp14(S);
Q(5)=Qcomp15(S);

% =================================================
function Q=Qcomp11(S)
% =================================================
   % Determine coefficients in the equation. 
    r=6*S(1)*S(4)^2*S(3)-4*S(2)^3*S(5)+6*S(3)^2*S(2)*S(4)+...
4*S(4)*S(2)^2*S(5)-6*S(3)^2*S(5)*S(2)-2*S(1)^2*S(4)*S(5)-...
6*S(1)*S(4)*S(3)^2-3*S(4)^2*S(3)^2+4*S(3)^3*S(5)+...
S(1)^2*S(5)^2+4*S(3)^3*S(1)+2*S(3)^3*S(2);

r=r+S(1)^2*S(4)^2+S(5)^2*S(2)^2-3*S(3)^2*S(2)^2-...
8*S(2)^2*S(4)^2+4*S(2)^3*S(4)-3*S(3)^4+2*S(3)^3*S(4)-...
4*S(4)^3*S(1)-2*S(5)^2*S(2)*S(1)-6*S(3)^2*S(1)*S(5);

r=r+6*S(3)*S(2)^2*S(5)+4*S(4)^2*S(1)*S(2)-...
6*S(5)*S(2)*S(4)*S(3)+6*S(4)*S(3)*S(1)*S(5)-...
6*S(1)*S(4)*S(3)*S(2)-2*S(1)*S(4)*S(5)*S(2)+...
6*S(3)*S(2)*S(1)*S(5)+4*S(4)^3*S(2);
Q=r;

% =================================================
function Q=Qcomp12(S)
% =================================================
r=8*S(1)*S(4)^2*S(3)-2*S(4)^2*S(3)*S(2)-2*S(6)*S(2)^2*S(5)+...
2*S(1)^2*S(4)*S(6)+6*S(6)*S(2)*S(3)^2+4*S(2)^3*S(6)-...
6*S(5)*S(2)*S(4)^2-4*S(1)*S(5)^2*S(3)-4*S(4)*S(2)^2*S(6)-...
6*S(1)*S(4)*S(3)^2-12*S(3)^2*S(5)*S(2)+4*S(5)^2*S(2)*S(3);

r=r+4*S(4)^2*S(1)*S(2)-2*S(1)^2*S(4)*S(5)-6*S(3)*S(2)^2*S(6)+...
12*S(4)*S(2)^2*S(5)+6*S(3)*S(2)^2*S(5)+6*S(3)^2*S(1)*S(6)-...
2*S(5)*S(2)*S(4)*S(3)+4*S(6)*S(2)*S(1)*S(5)-6*S(4)*S(3)*S(1)*S(6)+...
      2*S(4)*S(3)*S(1)*S(5)-10*S(1)*S(4)*S(5)*S(2);

r=r+2*S(1)*S(4)*S(6)*S(2)-6*S(3)*S(2)*S(1)*S(6)+6*S(3)*S(2)*S(1)*S(5)+...
6*S(1)*S(5)*S(4)^2-2*S(3)^4+6*S(6)*S(2)*S(4)*S(3)-6*S(4)^2*S(3)^2-...
4*S(3)^3*S(6)+4*S(3)^3*S(5)+4*S(4)^3*S(2)-4*S(2)^3*S(5)+...
2*S(1)^2*S(5)^2-6*S(4)^3*S(1)+2*S(3)^3*S(2)+6*S(3)^3*S(4)-...
2*S(5)^2*S(2)^2+2*S(4)^3*S(3)-4*S(2)^2*S(4)^2-2*S(1)^2*S(5)*S(6);

Q=r;
% =================================================
function Q=Qcomp13(S)
% =================================================

r=4*S(1)*S(5)^2*S(3)+6*S(5)^2*S(2)*S(4)+6*S(1)*S(5)*S(4)^2-...
  6*S(1)*S(5)^2*S(4)-6*S(3)^2*S(1)*S(5)-2*S(1)^2*S(5)*S(6)+...
6*S(1)*S(6)*S(4)^2-2*S(1)*S(6)^2*S(2)+6*S(4)*S(3)^2*S(6)+...
12*S(5)^2*S(2)*S(3)+6*S(3)^2*S(1)*S(6)-12*S(5)*S(2)*S(4)^2;

r=r-2*S(4)*S(3)^2*S(5)-14*S(4)^2*S(3)^2-3*S(3)^4-3*S(4)^4+...
2*S(1)*S(4)*S(6)*S(2)+2*S(1)*S(5)*S(6)*S(3)-2*S(5)*S(2)*S(6)*S(3)-...
6*S(4)^3*S(1)+2*S(6)*S(2)*S(1)*S(5)-14*S(4)*S(3)*S(1)*S(6)+...
2*S(4)*S(3)*S(1)*S(5)-2*S(1)*S(4)*S(5)*S(2);

r=r-8*S(3)^2*S(5)^2+4*S(4)*S(2)^2*S(6)+6*S(4)^2*S(3)*S(5)-...
     6*S(6)*S(2)*S(4)^2-12*S(3)^2*S(5)*S(2)-6*S(3)*S(2)^2*S(6)+...
6*S(3)^2*S(2)*S(4)+6*S(6)*S(2)*S(3)^2+6*S(3)*S(2)^2*S(5)+...
6*S(1)*S(4)^2*S(3)-2*S(4)^2*S(3)*S(2)+12*S(4)*S(2)^2*S(5);

r=r+S(1)^2*S(6)^2+12*S(3)^3*S(5)+12*S(4)^3*S(2)+S(1)^2*S(5)^2+...
6*S(3)^3*S(4)-9*S(5)^2*S(2)^2+6*S(4)^3*S(3)-8*S(2)^2*S(4)^2+...
S(6)^2*S(2)^2-6*S(3)^3*S(6)+2*S(6)*S(2)*S(4)*S(3)- ...
10*S(5)*S(2)*S(4)*S(3);

Q=r;
% =================================================
function Q=Qcomp14(S)
% =================================================


r=4*S(6)*S(3)^2*S(5)+2*S(5)*S(4)^3-4*S(3)^2*S(5)^2-4*S(4)*...
    S(2)^2*S(6)-2*S(5)^2*S(2)^2+6*S(1)*S(5)*S(4)^2-4*S(1)*S(5)^2*...
    S(3)-6*S(4)^2*S(3)*S(6)+4*S(3)^3*S(5)+2*S(3)^3*S(4);

r=r-4*S(4)^3*S(1)-...
2*S(5)*S(2)*S(4)*S(3)+4*S(6)*S(2)*S(1)*S(5)-6*S(4)*S(3)*S(1)*...
S(6)+6*S(4)*S(3)*S(1)*S(5)-6*S(1)*S(6)*S(5)*S(4)+2*S(1)*S(5)*S(6)*...
S(3)-10*S(5)*S(2)*S(6)*S(3)+6*S(5)*S(2)*S(6)*S(4)-2*S(4)^4+6*...
S(5)^2*S(2)*S(4);

r=r+ 2*S(6)*S(2)*S(4)*S(3)-6*S(3)^3*S(6)-6*S(4)^2*S(3)^2+4*...
     S(1)*S(5)^3+6*S(4)^3*S(3)+2*S(6)^2*S(2)^2-4*S(5)^3*S(2)+4*...
S(4)^3*S(2)+8*S(4)*S(3)^2*S(6)-2*S(4)*S(3)^2*S(5);

r=r-6*S(1)*S(5)^2*...
S(4)-2*S(6)^2*S(2)*S(3)-2*S(1)*S(6)^2*S(2)+6*S(1)*S(6)*S(4)^2-2*...
S(5)^2*S(2)*S(1)-6*S(3)^2*S(5)*S(2)+6*S(6)*S(2)*S(3)^2+2*S(1)*...
S(6)^2*S(3)+4*S(4)*S(2)^2*S(5)+12*S(5)^2*S(2)*S(3);
r=r-12*S(5)*S(2)*S(4)^2;

Q=r;
% =================================================
function Q=Qcomp15(S)
% =================================================

r=S(6)^2*S(3)^2-2*S(6)^2*S(2)*S(3)-6*S(6)*S(2)*S(4)^2+6*...
S(4)^2*S(3)*S(5)-2*S(6)*S(2)^2*S(5)-6*S(4)^2*S(3)*S(6)+6*S(4)*...
S(3)^2*S(6)+4*S(6)*S(3)^2*S(5)-6*S(5)*S(2)*S(4)^2+6*S(5)^2*S(2)*...
S(4)+4*S(5)^2*S(2)*S(3)+S(5)^2*S(2)^2-6*S(5)*S(2)*S(4)*S(3);

r=r+6*S(6)*S(2)*S(4)*S(3)-6*S(5)*S(4)*S(6)*S(3)+6*S(5)*S(2)*S(6)*S(4)-...
2*S(5)*S(2)*S(6)*S(3)-4*S(3)^3*S(6)+4*S(5)^3*S(3)-3*S(4)^2*S(3)^2-...
8*S(3)^2*S(5)^2+4*S(4)^3*S(6)+2*S(4)^3*S(3)-3*S(4)^4-4*S(5)^3*...
     S(2)+4*S(4)^3*S(2)-3*S(4)^2*S(5)^2+2*S(5)*S(4)^3+S(6)^2*...
S(2)^2+4*S(3)^3*S(5);

Q=r;

% =================================================
function R1=R1comp(S)
% =================================================

r=       -S(3)*S(2)^2-S(4)*S(2)^2-S(4)^2*S(1)-2*S(4)*S(3)*S(2)+2*...
         S(2)^2*S(5)-S(3)^3+3*S(3)^2*S(2)+S(4)*S(3)^2-2*S(3)*S(5)*...
         S(2)+S(4)^2*S(2)-S(1)*S(3)^2+2*S(3)*S(5)*S(1)+2*S(4)*S(2)*...
         S(1)-2*S(2)*S(5)*S(1);

R1(1)=r;

r=       -S(1)*S(3)^2-S(4)^2*S(1)-4*S(3)*S(5)*S(2)+3*S(2)^2*S(5)+...
         S(5)*S(3)^2-2*S(3)^3+4*S(3)^2*S(2)+2*S(4)*S(3)^2-2*S(4)*...
         S(2)^2-S(3)*S(4)^2+2*S(4)^2*S(2)-S(3)*S(2)^2-2*S(2)*S(5)*...
         S(1)+2*S(3)*S(5)*S(1)-2*S(4)*S(3)*S(2)+2*S(4)*S(2)*S(1);
R1(2)=r;

r=       -2*S(4)*S(3)*S(2)+3*S(3)^2*S(2)+2*S(3)*S(5)*S(1)-2*S(2)*S(5)*...
         S(1)+2*S(3)*S(4)^2-S(1)*S(5)^2-2*S(7)*S(3)*S(1)-6*S(3)*S(5)*...
         S(2)+3*S(5)*S(3)^2-2*S(3)*S(5)*S(4)-2*S(2)^2*S(7)+3*S(2)^2*...
         S(5)-S(4)^2*S(1)-2*S(3)^3-S(4)^3+2*S(4)^2*S(2)+2*S(2)*S(7)*...
         S(3)+S(2)*S(5)^2-S(4)*S(2)^2+2*S(5)*S(4)*S(1)+2*S(2)*S(7)*S(1);
R1(3)=r;

r=       2*S(2)^2*S(5)+2*S(7)*S(4)*S(1)+2*S(4)^2*S(2)-...
         6*S(3)*S(5)*S(4)+S(3)*S(5)^2+3*S(5)*S(4)^2-4*S(3)*S(5)*S(2)+3*S(5)*...
         S(3)^2+6*S(2)*S(7)*S(3)+2*S(4)*S(3)^2+2*S(3)*S(4)^2-2*S(4)^3-2*...
         S(7)*S(3)^2-2*S(2)*S(7)*S(4)-S(1)*S(5)^2-2*S(2)^2*S(7)-S(4)^2*...
         S(1)-2*S(4)*S(3)*S(2)+2*S(3)*S(5)*S(1)-2*S(7)*S(3)*S(1)-S(3)^3;
R1(4)=r;

r=      -2*S(2)*S(7)*S(4)-2*S(7)*S(3)^2+2*S(2)*S(7)*S(3)+S(2)*...
         S(5)^2-S(3)*S(4)^2-2*S(3)*S(5)*S(4)-2*S(4)*S(5)^2+3*S(5)*...
         S(4)^2+S(5)*S(3)^2-2*S(3)*S(5)*S(2)+S(4)^2*S(2)+2*S(7)*S(4)*...
         S(3)-S(4)^3+S(3)*S(5)^2+S(4)*S(3)^2;
R1(5)=r;

% =================================================
function Q=Qcomp2(S)
% =================================================

Q=zeros(5,5);

Q(1,1)= S(1)^2*S(4)^2+4*S(2)^3*S(4)-3*S(2)^2*S(3)^2+4*S(1)*S(3)^3-...
        6*S(2)*S(3)*S(1)*S(4);

Q(1,2)= -6*S(3)^2*S(1)*S(4)+4*S(1)*S(4)^2*S(2)-2*S(1)^2*S(4)*S(5)-...
        4*S(2)^3*S(5)+2*S(2)*S(3)^3+6*S(2)*S(3)*S(1)*S(5);

Q(1,3)= S(1)^2*S(5)^2-8*S(2)^2*S(4)^2+6*S(3)*S(4)^2*S(1)-2*S(1)*...
        S(4)*S(2)*S(5)+6*S(2)^2*S(3)*S(5)-3*S(3)^4+6*S(2)*S(3)^2*S(4)-6*...
        S(3)^2*S(1)*S(5);

Q(1,4)= -2*S(2)*S(5)^2*S(1)-6*S(3)^2*S(2)*S(5)+6*S(3)*S(4)*S(1)*S(5)+...
        2*S(3)^3*S(4)-4*S(1)*S(4)^3+4*S(2)^2*S(4)*S(5);

Q(1,5)= S(2)^2*S(5)^2-6*S(3)*S(4)*S(2)*S(5)-3*S(3)^2*S(4)^2+4*S(2)*...
        S(4)^3+4*S(3)^3*S(5);

%Q(2,1)=  0;

Q(2,2)= 8*S(3)*S(4)^2*S(1)+4*S(2)^3*S(6)-4*S(2)^2*S(4)^2-6*S(2)*S(3)*...
        S(6)*S(1)+6*S(2)^2*S(3)*S(5)-2*S(3)^4+2*S(1)^2*S(5)^2+2*S(1)^2*...
        S(4)*S(6)-10*S(1)*S(4)*S(2)*S(5);

Q(2,3)= 6*S(3)^3*S(4)+2*S(3)*S(4)*S(1)*S(5)-6*S(1)*S(4)^3+2*S(1)*S(4)*...
        S(6)*S(2)+6*S(3)^2*S(6)*S(1)-6*S(2)^2*S(3)*S(6)-2*S(4)^2*S(2)*...
        S(3)-12*S(3)^2*S(2)*S(5)-2*S(1)^2*S(5)*S(6)+12*S(2)^2*S(4)*S(5);

Q(2,4)= 6*S(4)^2*S(1)*S(5)+6*S(3)^2*S(6)*S(2)-6*S(3)*S(4)*S(6)*S(1)-...
        2*S(3)*S(4)*S(2)*S(5)-2*S(2)^2*S(5)^2-4*S(2)^2*S(4)*S(6)+4*S(2)*...
        S(5)*S(6)*S(1)+4*S(2)*S(4)^3-6*S(3)^2*S(4)^2+4*S(3)^3*S(5)-4*...
        S(1)*S(5)^2*S(3);

Q(2,5)= -6*S(4)^2*S(2)*S(5)-2*S(2)^2*S(5)*S(6)+6*S(3)*S(4)*S(6)*S(2)-...
        4*S(3)^3*S(6)+4*S(2)*S(5)^2*S(3)+2*S(4)^3*S(3);
%Q(3,1)=0;
%Q(3,2)=0;

Q(3,3)= 4*S(1)*S(5)^2*S(3)+2*S(2)*S(5)*S(6)*S(1)-9*S(2)^2*S(5)^2+12*...
        S(2)*S(4)^3+S(6)^2*S(1)^2-10*S(3)*S(4)*S(2)*S(5)+4*S(2)^2*S(4)*...
        S(6)-14*S(3)*S(4)*S(6)*S(1)+6*S(3)^2*S(6)*S(2)-14*S(3)^2*S(4)^2+...
        6*S(4)^2*S(1)*S(5)+12*S(3)^3*S(5);

Q(3,4)= 2*S(3)*S(4)*S(6)*S(2)+6*S(4)^2*S(6)*S(1)+12*S(2)*S(5)^2*S(3)-...
        6*S(1)*S(5)^2*S(4)-6*S(3)^3*S(6)-2*S(6)^2*S(2)*S(1)+2*S(1)*S(5)*...
        S(6)*S(3)+6*S(4)^3*S(3)-2*S(3)^2*S(5)*S(4)-12*S(4)^2*S(2)*S(5);

Q(3,5)= -8*S(3)^2*S(5)^2-6*S(4)^2*S(6)*S(2)+6*S(3)*S(4)^2*S(5)+6*...
        S(3)^2*S(4)*S(6)+S(6)^2*S(2)^2+6*S(2)*S(5)^2*S(4)-3*S(4)^4-2*...
        S(2)*S(5)*S(6)*S(3);

%Q(4,1:3)=zeros(1,3);

Q(4,4)= -10*S(2)*S(5)*S(6)*S(3)-6*S(6)*S(1)*S(5)*S(4)+6*S(2)*S(5)^2*...
        S(4)+2*S(6)^2*S(2)^2+8*S(3)^2*S(4)*S(6)+4*S(1)*S(5)^3-4*S(3)^2*...
        S(5)^2+2*S(6)^2*S(1)*S(3)-2*S(4)^4;

Q(4,5)= 2*S(4)^3*S(5)+4*S(3)^2*S(5)*S(6)+6*S(6)*S(2)*S(5)*S(4)-4*S(2)*...
        S(5)^3-2*S(6)^2*S(2)*S(3)-6*S(4)^2*S(6)*S(3);

%Q(5,1:4)=zeros(1,4);

Q(5,5)= 4*S(4)^3*S(6)-6*S(6)*S(3)*S(5)*S(4)-3*S(5)^2*S(4)^2+4*S(3)*...
        S(5)^3+S(6)^2*S(3)^2;

function z =f0min(Func,za, zb, options, R0, R1, R2, R3, R4, P, Q, D, sig)

Prob.R0  = R0;
Prob.R1  = R1;
Prob.R2  = R2;
Prob.R3  = R3;
Prob.R4  = R4;
Prob.P   = P;
Prob.Q   = Q;
Prob.D   = D;
Prob.sig = sig;

% z = fmin('exp_eq',za, zb, options, R0, R1, R2, R3, R4, P, Q, D, sig);                        
% Calling new TOMLAB version of fmin, called Tfmin
% z = Tfmin(Func, za, zb, options, Prob);


% xTol = 1E-7; Default is E-7, use that instead of Matlab fmin 1E-4
% z = tomsol(8, za, zb, Func, xTol, Prob);

z = tomsol(8, za, zb, Func, [], Prob);

% MODIFICATION LOG:
% 990622  hkh  Send empty array to f0min instead of calling foptions
% 011205  hkh  Qcomp2 should compute a 5x5 matrix, not a 25x1 vector
% 030116  hkh  Change comments to Tfmin instead of fmin
% 030118  hkh  Change name to Tfzero from dfzero
% 060814  med  FUNCS used for callbacks instead


% ========================================================================
%        exp_eq.m - local version with all parameters on the command line
% ========================================================================
%
% function y = exp_eq(z, R0, R1, R2, R3, R4, P, Q, D, sig);
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
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomopt.com
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.0.0$
% Written May 26, 1997.      Last modified Sep 30, 2000.
%

function y = exp_eq(z, R0, R1, R2, R3, R4, P, Q, D, sig)

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