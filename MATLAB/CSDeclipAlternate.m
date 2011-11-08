function [r,A] = CSDeclip(data)
%CSDECLIP - data is already clipped signal
%Naim Mansour

[rs cs]=size(data);
if(rs~=1)
    data=data';
end
[rs2 cs2]=size(data);
N=cs2;

X=[]; %support for reliable samples
E=[]; %support for missing samples
Mn=[];
Mp=[];
samples=[];
MaxS=max(data);
MinS=min(data);
for t=1:cs2  
    if(data(1,t)<max(data) && data(1,t)>min(data))
        X=[X t];
    elseif((t==1 || t==cs2)... 
            || (((length(find(MaxS)))==1) && data(1,t)==MaxS)...
            || (((length(find(MinS)))==1) && data(1,t)==MinS)...
            || (data(1,t-1)-data(1,t+1))==0)
        E=[E t];
        if data(1,t)==MaxS
            Mp=[Mp t];
        else
            Mn=[Mn t];
        end
    else
    end
end
samples=data';

%Remove all rows in the DCT base, according to the sampling matrix!
A=DCTBase(N,N,-1);
% A(E,:)=[]; %Important, removing the row, not making it zero!

if max(A)==0
    disp('Entire frame is clipped.')
    return;
end

%Remove all rows (NO, COLUMNS!!! -- in case of rows, do pinv(B)*B) in 
%the unit base, according to the error matrix
B=eye(length(samples),N);
I=eye(length(samples));
B(:,X)=[];

Re=I-B*pinv(B);
delta=eye(N);
for i=1:N
    delta(i,i)=1./norm(Re*A(:,i),2);
end
A=Re*A*delta;

%Solve the constrained L1 optimization (with lambda regularization)
% Delta_x=IRL1(A,samples,N,50,0.01,1e-4); %Development in progress
% Delta_x=SolveLasso(A,samples,N); --VERY SLOW, NOT THAT ACCURATE
% Delta_x=SolveBP(A,samples,N,50,0.01,1e-4); %Investigate parameter impact
Delta_x=SolveOMP(A,samples,N,50); %--FAST BUT INACCURATE

x=delta*Delta_x;

%HOW TO EXTRACT THE NOISE RECONSTRUCTION e?
%Seems a much worse approach.. ?
r=idct(x);

subplot(5,1,1);plot(data);
title('Clipped signal')
axis([0 N (min(data)-1) (max(data)+1)])
subplot(5,1,2);plot(r(1:N,1));
title('Reconstructed signal')
% subplot(5,1,3);plot(orig)
% title('Original signal')
% axis([0 N (min(orig)-1) (max(orig)+1)])
subplot(5,1,4);plot(samples)
title('Samples')
axis([0 N (min(samples)-1) (max(samples)+1)])
subplot(5,1,5);plot(abs(x))
title('Spectral representation of reconstruction')
end
