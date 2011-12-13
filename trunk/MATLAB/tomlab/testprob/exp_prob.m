% exp_prob.m
%
% Defines test problems for the fitting of positive sums of exponentials
% (nonlinear least squares fitting or separable nonlinear least squares)
%
% function [probList, Prob] = exp_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 7, 2008.

function [probList, Prob] = exp_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'Lanczos p=3 Artf'...
    ,'Lanczos p=3 Artf 4-dec'...
    ,'Lanczos p=3 Artf 2-dec'...
    ,'Lanczos p=4 Artf'...
    ,'EvansGV p=3 Artf 4-dec'...
    ,'EvansGV p=4 Artf 4-dec'...
    ,'EvansGV'...
    ,'EvansGV p=4 Artf'...
    ,'Generate artificial series'...
    ,'Boliden 1'...
    ,'Boliden 2'...
    ,'Boliden 3'...
    ,'Boliden 4'...
    ,'Astra 2'...
    ,'Astra 3'...
    ,'Steyn and Wyk 18points'...
    ,'Steyn and Wyk '...
    ,'Helax xp020 '...
    ,'Helax short xp020 '...
    ,'Helax non equidistant xp020 '...
    ,'Helax xp031 '...
    ,'Helax short xp031 '...
    ,'Helax non equidistant xp031 '...
    ,'Helax xp200 '...
    ,'Helax short xp200 '...
    ,'Helax non equidistant xp200 '...
    ,'Helax xp400 '...
    ,'Helax short xp400 '...
    ,'Helax non equidistant xp400 '...
    ,'Helax yp050 '...
    ,'Helax short yp050 '...
    ,'Helax non equidistant yp050 '...
    ,'Helax yp250 '...
    ,'Helax short yp250 '...
    ,'Helax non equidistant yp250 '...
    ,'Atexp nr1 '...
    ,'Atexp nr1\~ '...
    ,'Atexp nr2 '...
    ,'Atexp nr2\~ '...
    ,'Atcexp nr1 '...
    ,'Atcexp nr1\~ '...
    ,'Atcexp nr2 '...
    ,'Atcexp nr2\~ '...
    ,'TW_sim2'...
    ,'TW_sim5'...
    ,'TW4_real'...
    ,'TW3_sim2'...
    ,'CDFsim3'...
    ,'CDFdata3'...
    ,'CDFsim5'...
    ,'CDFdata5'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

SepAlg = []; x_0 = []; x_opt = [];
eType = 1;

if P==1
    Name='Lanczos p=3 Artf';
    t=(0:0.05:1.15)';
    y=0.0951*exp(-t)+0.8607*exp(-3*t)+1.5576*exp(-5*t);
    x_opt=[1 3 5];
elseif P==2
    Name='Lanczos p=3 Artf 4-dec';
    t=(0:0.05:1.15)';
    y=0.0951*exp(-t)+0.8607*exp(-3*t)+1.5576*exp(-5*t);
    x_opt=[1 3 5];
    y=0.0001*round(10000*y);
elseif P==3
    Name='Lanczos p=3 Artf 2-dec';
    t=(0:0.05:1.15)';
    y=[2.51 2.04 1.67 1.37 1.12 0.93 0.77 0.64 0.53 0.45 0.38 0.32 0.27 0.23 ...
        0.2 0.17 0.15 0.13 0.11 0.1 0.09 0.08 0.07 0.06]';
    x_opt=[1 3 5];
elseif P==4
    Name='Lanczos p=4 Artf';
    t=(1:0.2:32)';
    y=0.0951*exp(-t)+0.8607*exp(-3*t)+1.5576*exp(-5*t)+2*exp(-6*t);
    x_opt=[1 3 5 6];
elseif P==5
    Name='EvansGV p=3 Artf 4-dec';
    t=(0:1:17)';
    y=0.6*exp(-0.1*t)+0.3*exp(-0.01*t)+0.1*exp(-0.001*t);
    y=0.0001*round(10000*y);
    x_opt=[0.001 0.01 0.1];
elseif P==6
    Name='EvansGV p=4 Artf 4-dec';
    t=(1:0.2:32)';
    y=0.6*exp(-0.1*t)+0.3*exp(-0.01*t)+0.1*exp(-0.001*t)+0.2*exp(-0.0001*t);
    y=0.0001*round(10000*y);
    x_opt=[0.0001 0.001 0.01 0.1];
elseif P==7
    Name='Evans GV';
    t=[0 1 2 3 4 5 10 30 60 150 300 400 500 1000 1500 2000 3000 4000 5000 6000]';
    y=[1 0.9399 0.8853 0.8356 0.7904 0.7492 0.5921 0.3518 0.2655 0.1655 ...
        0.112 0.1016 0.09714 0.0905 0.08607 0.08187 0.07408 0.06703 ...
        0.06065 0.05488]';
elseif P==8
    Name='EvansGV p=4 Artf';
    t=(1:0.2:32)';
    y=0.6*exp(-0.1*t)+0.3*exp(-0.01*t)+0.1*exp(-0.001*t)+0.2*exp(-0.0001*t);
    x_opt=[0.0001 0.001 0.01 0.1];
elseif P==9
    % If not interactive mode, just define Exact Lanczos series instead.
    Name='Exact Lanczos p=3 Artf';
    t=(0:0.05:1.15)';
    y=0.0951*exp(-t)+0.8607*exp(-3*t)+1.5576*exp(-5*t);
elseif P==10
    Name='Boliden 1';
    t=[0 1 2.5 5 10 15 20 25 30 40]';
    y=[1 0.1543 0.075 0.05417 0.03528 0.02891 0.02472 0.02333 0.0225 ...
        0.02122]';
elseif P==11
    Name='Boliden 2';
    t=[0 2.25 5 10 15 20 25 30 40]';
    y=[1 0.1087 0.056 0.03397 0.02717 0.0216 0.01902 0.01821 0.01712]';
elseif P==12
    Name='Boliden 3';
    t=[0 1 3 5 10 15 20 25 30 40]';
    y=[1 0.1808 0.09973 0.07222 0.0486 0.03611 0.02861 0.02472 0.02111 ...
        0.01806]';
elseif P==13
    Name='Boliden 4';
    t=[0 1 3 5 10 15 20 25 30 40]';
    y=[1 0.1278 0.05361 0.03667 0.0208 0.01486 0.01139 0.01014 0.009167 ...
        0.008]';
elseif P==14
    Name='Astra 2';
    t = [0 0.033 0.083 0.17 0.25 0.333 0.5 0.75 1 1.5 2 2.5 3 4 5 6]' ;
    y = [29.655 45.716 40.323 28.868 24.643 20.05 14.81 10.869 7.674 ...
        4.136 2.776 1.587 1.135 0.511 0.294 0.2]';
elseif P==15
    Name='Astra 3';
    t = [0 0.03 0.08 0.17 0.25 0.33 0.5 0.75 1 1.5 2 2.5 3 4 5 6]';
    y = [22.583 37.34 33.981 25.534 24.013 17.405 13.934  9.039 6.933 ...
        4.176 2.947 1.886 1.37 0.691 0.412 0.24]';
elseif P==16
    Name='Steyn and Wyk 18points';
    t=(30:20:370)';   % Time in ms
    y=[18299 15428 13347 11466 10077 8729 7382 6708 5932 5352 4734 4271 ...
        3744 3485 3111 2950 2686 2476]';
    t=t/1000;     % Scale to seconds. Gives lambda*1000, of order 1
    y=y/10000;  % Scale function values. Avoid large alpha
elseif P==17
    Name='Steyn and Wyk ';
    t=(30:20:1170)';
    y=[18299 15428 13347 11466 10077 8729 7382 6708 5932 5352 4734 4271 ...
        3744 3485 3111 2950 2686 2476 2278 2107 1867 1768 1593 1504 1353 ...
        1389 1197 1080 1027 949 877 801 758 695 650 592 555 566 493 394 392 ...
        381 349 352 301 270 300 266 200 229 180 189 181 158 126 117 132 95]';
    t=t/1000;  % Scale to seconds. Gives lambda*1000, of order 1
    y=y/10000; % Scale function values. Avoid large alpha
    x_opt=[3.619756;11.548920;0.848085;1.507830];  % Found by tests, wType=0
elseif P>=18 & P <=35
    %Reads data (t, y) from Helax files.
    eType=2;
    Name=probList(P,:);
    i=1+floor((P-18)/3);
    j=mod(P-18,3);
    % ===================================================================
    % Helax Data Series ('xp020','xp031','xp200','xp400','yp050','yp250')
    % ===================================================================
    if i==1
        load xp020 A
    elseif i==2
        load xp031 A
    elseif i==3
        load xp200 A
    elseif i==4
        load xp400 A
    elseif i==5
        load yp050 A
    elseif i==6
        load yp250 A
    end
    %A=feval(HFile);  % Define matrix A with Helax data
    n=size(A,1);
    if j==0
        t=A(:,1);
        y=A(:,2);
    elseif j==1   % Define short series
        t=A(1:min(20,n),1);
        y=A(1:min(20,n),2);
    else          % Define non equidistant series
        ix=ones(n,1);
        ix(3:3:n)=0;
        t=A(ix>0,1);
        y=A(ix>0,2);
    end
elseif P==36
    Name='Atexp nr1 ';
    t=(0:0.1:1.7)';
    y=t.*exp(-0.1*t)+2*t.*exp(-0.2*t);
    x_opt=[0.1 0.2];
    eType=3;
elseif P==37
    Name='Atexp nr1\~ ';
    t=(0:0.1:1.7)';
    y=t.*exp(-0.1*t)+2*t.*exp(-0.2*t);
    y=0.0001*round(10000*y);
    eType=3;
elseif P==38
    Name='Atexp nr2 ';
    t=(0:0.1:1.7)';
    y=1.1*t.*exp(-0.08*t)+2.5*t.*exp(-0.22*t);
    x_opt=[0.08 0.22];
    eType=3;
elseif P==39
    Name='Atexp nr2\~ ';
    t=(0:0.1:1.7)';
    y=1.1*t.*exp(-0.08*t)+2.5*t.*exp(-0.22*t);
    y=0.0001*round(10000*y);
    eType=3;
elseif P==40
    Name='Atcexp nr1 ';
    t=(0:0.1:1.7)';
    y=(t+0.5).*exp(-0.1*t)+(2*t+1).*exp(-0.2*t);
    x_opt=[0.1 0.2];
    eType=4;
elseif P==41
    Name='Atcexp nr1\~ ';
    t=(0:0.1:1.7)';
    y=(t+0.5).*exp(-0.1*t)+(2*t+1).*exp(-0.2*t);
    x_opt=[0.1 0.2];
    y=0.0001*round(10000*y);
    eType=4;
elseif P==42
    Name='Atcexp nr2 ';
    t=(0:0.1:1.7)';
    y=(1.1*t+0.6).*exp(-0.08*t)+(2.5*t+0.2).*exp(-0.22*t);
    x_opt=[0.08 0.22];
    eType=4;
elseif P==43
    Name='Atcexp nr2\~ ';
    t=(0:0.1:1.7)';
    y=(1.1*t+0.6).*exp(-0.08*t)+(2.5*t+0.2).*exp(-0.22*t);
    x_opt=[0.08 0.22];
    y=0.0001*round(10000*y);
    eType=4;
elseif P==44
    % Simulation for type 2 problem
    Name='TW_sim2';
    rand('state',0);
    [t,y,yModel]=WeibSim(0.7,[	1 0.5],zeros(2,1),200);
    y=yModel;
    dig=0; % No not round now.
    if dig>0, t=round(dig*t)/dig; end
    eType = 2;
elseif P==45
    % Simulation for type 5 problem
    Name='TW_sim5';
    rand('state',0);
    %HKH Better formulated problem, no negative values:
    [t,y]=WeibSim(0.7,[0.5 0.001],[0.5 0.1],200);
    dig=10000;   % Round to five decimals in t
    if dig>0, t=round(dig*t)/dig; end
    eType = 1;
elseif P==46
    % Real data Smax for yet to be created eType=5 problem.
    Name='TW4_real';
    t=[0.56 0.67 0.68 0.68 0.71 0.72 0.73 0.74 0.77 0.78 0.81 0.82 0.84...
        0.84 0.86 0.88 0.89 0.89 0.97 0.99 ...
        0.99 1.01 1.02 1.02 1.03 1.04 1.04 1.04 1.05 1.05 1.06 1.07 1.09...
        1.13 1.13 1.13 1.17 1.21 1.25 1.28 ...
        1.35 1.36 1.36 1.40 1.41 1.42 1.43 1.43 1.44 1.44 1.47 1.48 1.57 ...
        1.63 1.71 1.72 1.73 2.15 2.28 2.48]';
    %  Real data with [b1 b2 c1 c2 a1 a2] starting values ~
    %                 [2.5044 1.6252 0.8116 0.5497 0.5 0.5]
    y=(1:60)/61;y=y';
    eType = 1;
elseif P==47
    % Simulation for real type 2 problem
    Name='TW3_sim2';
    t=[0.009 0.015 0.018 0.039 0.059 0.149 0.193 0.221 0.2215 0.226 ...
        0.262 0.317 0.420 0.435 0.452 0.520 0.528 0.588 0.5885 0.609 ...
        0.627 0.665 0.726 0.925 0.933 0.954 0.955 1.084 1.120 1.340 ...
        1.373 1.397 1.435 1.489 1.569 1.677 1.722 2.219 2.230 2.241 ...
        2.286 2.472 2.487 2.548 2.685 2.740 2.998 3.565 3.641 3.744]';
    y=(1:50)/51;y=y';
    eType = 2;
elseif P==48
    % eType=5 problem , constrained exponential fit with 2 exponentials
    % yModel=F(x) [i.e. perfect fit] OR
    % with y=empirical distribution F(i/(N+1))
    % F(x)=a1*(1-exp(-b1*t))+(1-a1)*(1-exp(-b2*t))
    Name='CDFsim3';
    rand('state',0);
    %[t,y,yModel]=WeibSim(0.7,[	1 0.5],[0 0],200);
    % changed zeros(2,1) to [0 0]
    [t,y,yModel]=WeibSim(0.7,[3.0 0.3],[0 0],200);
    y=yModel;
    dig=0; % No not round now.
    if dig>0, t=round(dig*t)/dig; end
    eType = 1;
elseif P==49
    % eType=5 problem , constrained exponential fit with 2 exponentials
    % F(t)=a1*(1-exp(-b1*t))+(1-a1)*(1-exp(-b2*t))
    Name='CDFdata3';
    x=[0.56 0.67 0.68  0.685 0.71  0.72  0.73  0.74  0.77 0.78 ...
        0.81 0.82 0.84  0.845 0.86  0.88  0.89  0.90  0.97 0.99 ...
        0.99 1.01 1.02  1.025 1.03  1.04  1.045 1.048 1.05 1.055 ...
        1.06 1.07 1.09  1.13  1.135 1.138 1.17  1.21  1.25 1.28 ...
        1.35 1.36 1.365 1.40  1.41  1.42  1.43  1.435 1.44 1.445 ...
        1.47 1.48 1.57  1.63  1.71  1.72  1.73  2.15  2.28 2.48 ]';
    t=x-min(x);
    y=(1:length(t))/(length(t)+1); y=y';
    dig=0; % No not round now.
    if dig>0, t=round(dig*t)/dig; end
    eType = 1;
elseif P==50
    % eType=5 problem , constrained exponential fit with 2 exponentials
    % yModel=F(x) [i.e. perfect fit]    OR       with y=empirical
    %         distribution F(i/(N+1))
    % F(x)=a1*(1-exp(-b1*(t-c1)))+(1-a1)*(1-exp(-b2*(t-c2)))
    Name='CDFsim5';
    rand('state',0);
    [t,y]=WeibSim(0.7,[1 1.5],[0.5 0.3],200);
    dig=10000;   % Round to five decimals in t
    if dig>0, t=round(dig*t)/dig; end
    eType = 1;
elseif P==51
    % eType=5 problem , constrained exponential fit with 2 exponentials
    Name='CDFdata5';
    x=[0.56 0.67 0.68  0.685 0.71  0.72  0.73  0.74  0.77 0.78 ...
        0.81 0.82 0.84  0.845 0.86  0.88  0.89  0.90  0.97 0.99 ...
        0.99 1.01 1.02  1.025 1.03  1.04  1.045 1.048 1.05 1.055 ...
        1.06 1.07 1.09  1.13  1.135 1.138 1.17  1.21  1.25 1.28 ...
        1.35 1.36 1.365 1.40  1.41  1.42  1.43  1.435 1.44 1.445 ...
        1.47 1.48 1.57  1.63  1.71  1.72  1.73  2.15  2.28 2.48 ]';
    t=x;
    y= (1:length(t))/(length(t)+1); y=y';
    eType = 1;
else
    error('exp_prob: Illegal problem number');
end

p = 2;
weightType=3;
Prob = expAssign(p, Name, t, y, weightType, eType, SepAlg, x_0);
Prob.x_opt = x_opt;
Prob.P = P;

% ===================================================================
% Internal functions exp_ArtP, WeibInit, WeibIni2 and WeibSim
% ===================================================================

% ===================================================================
%		exp_ArtP.m
% ===================================================================
%
% function [Z, prob] = exp_ArtP (p,lambda,alpha,t,m,distP,PriLev);
%
% Generate m series of artificial data for exponential approximation problem
%
% INPUT:  p      # of exponential terms. Default p=2.
%         lambda exponential intensities. Default [1,2,...,p].
%         alpha  weights. Default 100*ones(p,1);
%         t      A N-vector with time points. Default [0.1:0.1:3]
%         m      # of series to generate. Default 1.
%         distP  Parameters determining the generated random data
%            (1) Random variable distribution:
%                1 = Uniformly distributed random numbers
%                2 = Random numbers from normal distribution
%                3 = Random numbers from exponential distribution
%            (2) Standard deviation to multiply the random value with.
%            (3) Error type:
%                1 = Additive error on log y
%                2 = Relative error on log y
%                3 = Additive error on y
%         PriLev 0=Silent. >1 Plot each time series
%
% OUTPUT: Z      Z is N x m, where each col is an artifical time series
%         prob   Name of generated problem. String.
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Written Mar 11, 1997.  Last modified June 21, 1999.
%

function [Z, prob] = exp_ArtP(p,lambda,alpha,t,m,distP,PriLev)

if nargin < 7
    PriLev=[];
    if nargin < 6
        distP=[2 0.2 1];
        if nargin < 5
            m=[];
            if nargin < 4
                t=[];
                if nargin < 3
                    alpha=[];
                    if nargin < 2
                        lambda=[];
                        if nargin < 1
                            p=[];
                        end
                    end
                end
            end
        end
    end
end

if isempty(p)==0,      p=2;                 end
if isempty(lambda)==0, lambda=(1:p)';       end
if isempty(alpha)==0,  alpha=100*ones(p,1); end
if isempty(t)==0,      t=(0.1:0.01:3)';     end
if isempty(m)==0,      m=1;                 end
if isempty(distP)==0,  distP=[2 0.2 1];     end
if isempty(PriLev)==0, PriLev=2;            end
if length(distP)<3
    dp=[2 0.2 1];
    dp(1:length(distP))=distP;
    distP=dp;
end

t=t(:);
lambda=lambda(:);
alpha=alpha(:);
N=length(t);
Z=zeros(N,m);
El=exp(-t*lambda');
y=El*alpha;

if PriLev > 2
    clf;
    title('Exponential time series'); hold on;
    xlabel('t');
    plot(t,y);
    pause
end

% Generate N*m random values r
if distP(1)==1
    r=rand(N, m);
elseif distP(1)==2
    r=randn(N, m);
elseif distP(1)==3
    r=rand(N, m);
    r=0.5*log(1./(1-r));    %exp=f^(-1)(unif). See Blom G. page 217.
end

% Compute artificial series
if distP(3)==1
    Z=exp(log(y)*ones(1,m)+r*distP(2));
elseif distP(3)==2
    Z=exp(log(y)*ones(1,m).*(1+r*distP(2)));
elseif distP(3)==3
    Z=y*ones(1,m)+r*distP(2);
elseif distP(3)==4
    Z=y*ones(1,m).*(1+r*distP(2));
end

% Generate name of problem
prob=['Art Prob p=' num2str(p') ' lambda='];
for i=1:p
    prob=[prob ' ' num2str(lambda(i))];
end
prob=[prob ' alpha='];
for i=1:p
    prob=[prob ' ' num2str(alpha(i))];
end

if PriLev > 0
    disp(prob);
end

if PriLev > 1
    for i=1:m
        clf;
        title(['Artificial and original time series #' num2str(i)]);
        hold on;
        xlabel('t');
        %plot(t,[Z(:,i),y]);
        plot(t,y,'r');
        plot(t,Z(:,i),'x');
        if PriLev > 4, pause, end
    end
end

% ===================================================================
% function WeibInit(x);
% ===================================================================
% Calculates the parameters for a shifted exponential on
% partial data set (x=data vector);
% May be used as initial values for an parameter optimization
%
% Model:
% y=a1(1-exp(-b1(t-c1)))+(1-a1)(1-exp(-b2(x-c2)));
%
% It simply fits the lower portion of data to t=x(i), y=F(x(i)) where
% F(x)=1-exp(-b(x-c));
% Since using half the data, a=0.5.
% After fitting lower half of data, it fits the upper half in a similiar manner
% Program gives values of b,c for each exponential.
%
% Written   by Todd Walton       Feb 18, 1999.
% Rewritten by Kenneth Holmstrom Feb 21, 1999.

function [a,b,c] = WeibInit(x)

xsort=sort(x);
N=length(x);

ivec=1:N;ivec=ivec';
Femp=ivec/(N+1);

% estimate parameters for 1 exponental on entire set of data
c_1=min(x);
beta_ML=sum(x-c_1)/(N-1);
x0_ML=c_1-beta_ML/N;

% choose portion of data to fit
xt=xsort(round(N/2));

a = [round(N/2)/N, 1-round(N/2)/N];

%    fit lower half first
isub=find(xsort<xt); xsub=xsort(isub); Fe_sub=Femp(isub);
%   initial guess of para
para0=[beta_ML; x0_ML];
%
Prob = conAssign('funexpm2',[],[],[],[],[],'weibinit',para0);
Prob.user.xSub=xsub;
Prob.user.FeSub=Fe_sub;
Result = tomRun('ucSolve',Prob,0);
para = Result.x_k;

% NOT USE fmins
%para=fmins('funm2exp',para0,[],[],xsub,Fe_sub);

betalow=para(1);x0low=para(2);
c_low=x0low;   % c for exponential 1
b_low=1/betalow; % b for exponential 1

% a for exponential 1 =0.5

%    fit upper half next
isub=find(xsort>xt); xsub=xsort(isub); Fe_sub=Femp(isub);
%   initial guess of para
para0=[beta_ML; x0_ML];

Prob.user.xSub=xsub;
Prob.user.FeSub=Fe_sub;
Prob.x_0 = para0;
Result = tomRun('ucSolve',Prob,0);
para = Result.x_k;

% NOT USE fmins
%para=fmins('funm2exp',para0,[],[],xsub,Fe_sub);

betahigh=para(1);x0high=para(2);
c_high=x0high;  % c for exponential 2
b_high=1/betahigh; % b for exponential 2
% a for exponential 2 = 0.5

b=[b_low b_high];
c=[c_low c_high];

% =====================================================================
% function [x_0]=WeibIni2(t);
% =====================================================================
% modified from WeibInit.m
% Calculates the initial guess of parameters for eType==5 problem with
%   c1,c2=0;
% May be used as initial values for an parameter optimization
%
% Model:
% y=a1(1-exp(-b1*t))+(1-a1)(1-exp(-b2*t)); where t=sorted vector
%                                                  =sort(x-min(x));
%
% It fits the two halves of data with separate exponential functions;
% Since using half the data, a=0.5.
%
%     x_0= initial estimate of parameters=[b1 b2 0 0 0.5 0.5]
%                                        =[b1 b2 c1 c2 a1 (1-a1)];
%
% Written   by Todd Walton       Feb 18, 1999.
% Rewritten by Kenneth Holmstrom Feb 21, 1999.
% Rewritten by Todd Walton       Feb 22, 1999.

function [x_0] = WeibIni2(t)
%
N=length(t);
t=sort(t); % a precaution
% Recall t=sort(x-xmin);
ivec=1:N;ivec=ivec';
Femp=ivec/(N+1);

% estimate parameters for 1 exponental on entire set of data
tmid=median(t);

% Lets do it via a linear approach instead;
disp('Linear approach to calculate initial parameters');
Y=-log(ones(N,1)-Femp);
%    fit lower half first
isub=find(t<=tmid); tsub=t(isub); Y_sub=Y(isub);A=tsub;
b_low=pinv(A'*A)*(A'*Y_sub);
%    fit upper half now;
isub=find(t > tmid); tsub=t(isub); Y_sub=Y(isub);A=tsub;
b_high=pinv(A'*A)*(A'*Y_sub);
% set up parameter vector- Recall this is for t=xsort-min(x);
x_0=[b_low b_high 0 0 0.5 0.5];  % input x_0 initial vector to be
%                    included in 'exp_prob',P==44.

% ===================================================================
% function [t, y, yModel]=WeibSim(a, b, c, N);
% ===================================================================
%
% Generate eType = 2 or eType = 5 type of simulated data
%
% f(t)=sum(i=1:p) a(i)*(1-exp(-b(i)*(x-c(i))))
%
% e.g for two terms
% f(t)=a1*(1-exp(-b1*(x-c1)))+a2*(1-exp(-b2*(x-c2)))
%
% Number of terms p is computed as p = length(a)+1
%
% a   Only give the first p-1 elements of a, a(p)=1-sum(i=1:p-1) a(i)
% b,c Must have length(p)
%
%
% N    length of data series (default 50)
%
% Written   by Todd Walton       Feb  5, 1999.
% Rewritten by Kenneth Holmstrom Feb 16, 1999.

function [t, y, yModel]=WeibSim(a, b, c, N)

if nargin < 4
    N=50;     %length of series desired
    if nargin < 3
        c=[0.5 1.5];
        if nargin < 2
            b=[1 1];
            if nargin < 1
                a=0.7;
            end
        end
    end
end

a=a(:); b=b(:); c=c(:);

a=[a;1-sum(a)];
p=length(a);
x=[];
n=round(a*N);
n(p)=N-sum(n(1:p-1));  % Avoid errors in length due to rounding

for i=1:p
    x=[x;c(i)-log(1-rand(n(i),1))/b(i)];
end

t=sort(x);
y=((1:N)/(N+1))'; % y for Weibull plotting position assumption

% The model
E = exp((t*ones(1,p)-ones(N,1)*c')*diag(-b));
yModel=(1-E)*a;
ix=find(yModel < 0);

if ~isempty(ix)
    fprintf('WeibSim: Warning. yModel < 0 for %d points\n',length(ix));
    xprinte(yModel(ix),'<0:');
end

% ===================================================================
function Prob = WeibProb(Prob,P,t)
% ===================================================================
p=Prob.ExpFit.p;
% Avoid too many iterations now when we have problem with this model
Prob.optParam.MaxIter=100;

if P < 48 | P > 49
    %KH put in to calc initial estimate via modified TLW kh_init.m
    % NOT USED NOW
    %[a,b,c]=WeibInit(t);
else
    [x_0]=WeibIni2(t); % TLW file for x_0=[b1 b2 0 0 0.5 0.5];
end

if P==44
    Prob.x_0=[1 0.5 0 0 0.7 0.3]';
    % Start with the next value if the run yModel instead of y
    %Prob.x_0=[1.2 0.6 0.1 0.2 0.5 0.5]';
elseif P==45
    %Starting values for:
    %[t,y,yModel]=WeibSim(0.7,[0.5 0.001],[0.5 0.1],200);
    % Start at optimum
    %Prob.x_0=[0.5 0.001 0.5 0.1 0.7 0.3]';
    % Start far away, problem with c to converge
    Prob.x_0=[0.7 0.010 0.4 0.2 0.5 0.5]';

    %Starting values for:
    %[t,y,yModel]=WeibSim(0.7,[1 1.5],[0.5 0.3],200);
    %[t,y,yModel]=WeibSim(0.7,[0.5 0.001],[0.5 0.1],200);
    %Prob.x_0=[1 1.5 0.5 1.5 0.7 0.3]';
elseif P==46
    % Running the starting value below clsSolve stops to early
    %Prob.x_0=[1 1 0.5 1.5 0.7 0.3]';
    % With this starting value it converges fast
    Prob.x_0=[1 0.5 0 0 0.7 0.3]';
    % It takes a while to converge, and it gets stuck
    %Prob.x_0=[1:p,zeros(1,p),ones(1,p)/p]';
    % The starting values supplied by Todd Walton
    Prob.x_0=[2.5044 1.6252 0.8116 0.5497 0.5 0.5]';
elseif P==48
    Prob.x_0=x_0'; % obtain values from WeibIni2.m
elseif P==49
    Prob.x_0=x_0'; % obtain values from WeibIni2.m
elseif P==50
    %Starting values for:
    %[t,y,yModel]=WeibSim(0.7,[1 1.5],[0.5 0.3],200);
    Prob.x_0=[1 1.5 0.5 0.3 0.7 0.3]';
elseif P==51
    % Running the starting value below clsSolve stops to early
    %Prob.x_0=[1 1 0.5 1.5 0.7 0.3]';
    % With this starting value it converges fast
    Prob.x_0=[1 0.5 0 0 0.7 0.3]';
    % It takes a while to converge, and it gets stuck
    %Prob.x_0=[1:p,zeros(1,p),ones(1,p)/p]';
    %Prob.x_0=[2.5044 1.6252 0.8116 0.5497 0.5 0.5]';
    tomRun('nlssol', Prob);
    % start values supplied by TW
else
    Prob.x_0=[1:p,zeros(1,p),ones(1,p)/p]';
end
if length(Prob.x_0) ~= 3*p
    x0  = Prob.x_0(1:2);
    my0 = Prob.x_0(3:4);
    a0  = Prob.x_0(5:6);
    if p == 1
        Prob.x_0 = [x0(1);my0(1);a0(1)];
    else
        Prob.x_0 = [x0; x0(2)+2*(1:p-2)'; my0; zeros(p-2,1); ...
            a0; 1E-4*ones(p-2,1)];
    end
end
Prob.x_L=[zeros(1,p) zeros(1,p) zeros(1,p)]';
Prob.x_min=[zeros(1,p) zeros(1,p) zeros(1,p)]';
Prob.x_U=[50*ones(1,p) min(t)*ones(1,p) 1.2*ones(1,p)]';
Prob.x_0=min(Prob.x_U,Prob.x_0);
Prob.x_max=[10*ones(1,p) min(t)*ones(1,p) 1.0*ones(1,p)]';

if Prob.LS.SepAlg
    % lambda ,beta parameters
    Prob.x_0=Prob.x_0(1:2*p);
    Prob.x_L=Prob.x_L(1:2*p);
    Prob.x_U=Prob.x_U(1:2*p);
    Prob.x_min=Prob.x_L;
    Prob.x_max=Prob.x_U;
end
if ~Prob.LS.SepAlg
    % The sum of the weights should be 1
    Prob.A=[zeros(1,2*p) ones(1,p)];
    Prob.mLin = size(Prob.A,1);
    Prob.b_L=1;
    Prob.b_U=1;
end

% MODIFICATION LOG
%
% 980825  hkh  Check if P is a field, otherwise initialize.
%
% 980826  hkh  Defining probType before call to ProbVarDef.
% 980922  hkh  Change name f_min to f_Low
% 981022  mbk  c_L=[]; c_U=[]; for all problems.
% 981022  hkh  Changed to use tomFiles for name definitions
% 981022  hkh  Set Name,t,y=[] if P is illegal and return
% 981022  hkh  Made safe and general way to read Helax data files.
% 981023  hkh  Use ls_H for Hessian
% 981027  hkh  Check which P in Prob.P, not just on if nonempty
% 981028  hkh  Define linear constraints in new routine ExpLinCon instead
%              Add artificial problems for eType=3,4
% 981102  hkh  Put exact Lanczos as problem 9, if ask~=1, to avoid questions
%              in batch runs
% 981105  hkh  Safeguard test on Prob.P and Prob.probFile
% 990204  hkh  Add from P=44 from Todd Walton
% 990216  hkh  Expand routine to handle eType=5; distribution estimation
%              Add eType=5 problems.
% 990322  hkh  Add more problems from Todd Walton
% 990623  hkh  Use clsVarDef and clsProbSet instead of ProbVarDef and ProbSet
% 020222  hkh  Change use of ExpLinCon and use it externally defined
% 031120  med  Add ; for CDFdata5 at many places, avoid output
% 040126  hkh  Added definition of field mLin
% 041117  med  xxx_prob removed and code added
% 080607  hkh  Use tomRun not tomSolve in Todd Walton example
