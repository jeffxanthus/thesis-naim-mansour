function [ threshold ] = maskingThreshold( input )
%calculates masking threshold of input audio signal
%   Steven De Hertogh
% hoi


% Reading inputsignal
[in, Fs, nbits] = wavread(input);

%FFT length
N = 512;

% Normalization
sig = [];
for i = 1:length(in)
sig(i) = in(i,1)/(N*2^(nbits-1));
end
%sig = sig(10000:50000);

figure();
subplot(3,1,1);plot(in(:,1));
%subplot(3,1,2);plot(in(:,2));
subplot(3,1,3);plot(sig);

% PSD
P=[];
z=zeros(1,N/2);
PN = 90;
for k = 1:N/2
    w = [];
   for n = 1:N
       w(n) = 1/2*(1-cos(2*pi*n/N));
       z(k) = z(k) + w(n)*in(n,1)*exp(-j*2*pi*k*n/N);
   end
   P(k) = PN + 10*log10(abs(z(k))^2);
end
figure();
plot(P);

% Tonal set
St=zeros(1,N/2);

for k = 3:62 %(0,17-5,5 kHz)
    if (P(k) > P(k+1)) && (P(k) > P(k-1)) && (P(k) > (P(k+2)+7)) && (P(k) > (P(k-2)+7))
         St(k) = P(k);
    end
end
for k = 63:126 %(5,5-11 kHz)
    if (P(k) > P(k+1)) && (P(k) > P(k-1)) && (P(k) > (P(k+2)+7)) && (P(k) > (P(k-2)+7)) && (P(k) > (P(k+3)+7)) && (P(k) > (P(k-3)+7))
         St(k) = P(k);
    end
end
for k = 127:250 %(11-20 kHz) 
    if (P(k) > P(k+1)) && (P(k) > P(k-1)) && (P(k) > (P(k+2)+7)) && (P(k) > (P(k-2)+7)) && (P(k) > (P(k+6)+7)) && (P(k) > (P(k-6)+7))
         St(k) = P(k);
    end
end

i = 1;
a= [];
b= [];
for k = 1:256
   if St(k) ~= 0
       a(i) = St(k);
       b(i) = k;
       i = i+1;
   end
end
hold on
plot(P);
plot(b,a, 'o');
%plot(St, 'o');
hold off

% Tonal and Noise Maskers
Ptm = zeros(1,N/2);

for k = 2:(N/2-1)
   w = 0;
   if St(k) ~= 0
    for c = -1:1
       w = w + 10^(0.1*St(k+c));
    end
   end
   Ptm(k) = 10*log10(w);
end

Pn = zeros(1,N/2);
for k = 3:62 %(0,17-5,5 kHz)
    if (P(k) ~= St(k)) && (P(k) ~= St(k+1)) && (P(k) ~= St(k-1)) && (P(k) ~= St(k+2)) && (P(k) ~= St(k-2))
         Pn(k) = P(k);
    end
end
for k = 63:126 %(5,5-11 kHz)
    if (P(k) ~= St(k)) && (P(k) ~= St(k+1)) && (P(k) ~= St(k-1)) && (P(k) ~= (St(k+2)+7)) && (P(k) ~= St(k-2)) && (P(k) ~= St(k+3)) && (P(k) ~= St(k-3))
         Pn(k) = P(k);
    end
end
for k = 127:250 %(11-20 kHz) 
    if (P(k) ~= St(k)) && (P(k) ~= St(k+1)) && (P(k) ~= St(k-1)) && (P(k) ~= (St(k+2)+7)) && (P(k) ~= St(k-2)) && (P(k) ~= St(k+6)) && (P(k) ~= St(k-6))
         Pn(k) = P(k);
    end
end

Pnm = zeros(1,N/2);
i = 1;
k = 1;
while k < N/2
    kstart = k;
    while ((ceil(FreqToBark(k*86)) == i) && (k < (N/2)))
        k = k+1;
    end
    kstop = k;
    i = i+1;
    kstreep = 1;
    for z = kstart:kstop
       kstreep = kstreep * z;
    end
    kstreep = kstreep^(1/(kstop-kstart+1));
    kstreep = floor(kstreep);
    w = 0;
    nottonal = 1;
    for n = kstart:kstop
        if Pn(n) == 0
            nottonal = 0;
        else
            w = nottonal*(w + 10^(0.1*Pn(n)));
        end
    end
    Pnm(kstreep) = 10*log10(w);
end


i = 1;
n = 1;
a= [];
b= [];
c= [];
d= [];
for k = 1:256
   if Ptm(k) ~= 0
       a(i) = Ptm(k);
       b(i) = k;
       i = i+1;
   end
   if Pnm(k) ~= 0
       c(n) = Pnm(k);
       d(n) = k;
       n = n+1;
   end
       
end
figure();
hold on
plot(P);
plot(b,a, 'o');
plot(d,c, 'x');
hold off





% Absolute threshold
Tq = [];
for f = 1:22000
    Tq(f) = 3.64*(f/1000)^(-0.8) - 65*exp(-0.6*(f/1000-3.3)^2) + 10^(-3) * (f/1000)^4;
end

tq=(downsample(Tq, 86));

figure();
hold on
plot(P);
plot(b,a, 'o');
plot(d,c, 'x');
semilogx(tq);
hold off


% Decimation of maskers
for k = 1:N/2
    if Ptm(k) < tq(k)
        Ptm(k) = 0;
    end
    if Pnm(k) < tq(k)
        Pnm(k) = 0;
    end
end



i = 1;
n = 1;
a= [];
b= [];
c= [];
d= [];
for k = 1:256
   if Ptm(k) ~= 0
       a(i) = Ptm(k);
       b(i) = k;
       i = i+1;
   end
   if Pnm(k) ~= 0
       c(n) = Pnm(k);
       d(n) = k;
       n = n+1;
   end
end
figure();
hold on
plot(P);
plot(b,a, 'o');
plot(d,c, 'x');
semilogx(tq);
hold off



% Sliding 0.5 bark window
%tone
ks = [];
b = [];
i=1;
for k = 1:N/2
   if Ptm(k)>0
       ks(i)=k;
       b(i) =  FreqToBark(k*86);
       i = i+1;
   end    
end


max = 0;
maxki = 1;
for i = 1:length(ks)
     if (b(i) - b(maxki)) < 0.5
         if (Ptm(ks(i)) > max)
             max = Ptm(ks(i));
             maxki = i;
         else
             Ptm(ks(i)) = 0;
         end
     else
     max = 0;
     maxki = i;
     end
end

%noise
ks = [];
b = [];
i=1;
for k = 1:N/2
   if Pnm(k)>0
       ks(i)=k;
       b(i) =  FreqToBark(k*86);
       i = i+1;
   end    
end


max = 0;
maxki = 1;
for i = 1:length(ks)
     if (b(i) - b(maxki)) < 0.5
         if (Pnm(ks(i)) > max)
             max = Pnm(ks(i));
             maxki = i;
         else
             Pnm(ks(i)) = 0;
         end
     else
     max = 0;
     maxki = i;
     end
end





i = 1;
n = 1;
a= [];
b= [];
c= [];
d= [];
for k = 1:256
   if Ptm(k) ~= 0
       a(i) = Ptm(k);
       b(i) = k;
       i = i+1;
   end
   if Pnm(k) ~= 0
       c(n) = Pnm(k);
       d(n) = k;
       n = n+1;
   end
end
figure();
hold on
plot(P);
plot(b,a, 'o');
plot(d,c, 'x');
semilogx(tq);
hold off


%Decimation and reorganization
Ptm2 = zeros(1,length(Ptm));
Pnm2 = zeros(1,length(Pnm));
for k = 1:length(Ptm)
    if  (k >= 1) && (k <= 48)
        i = k;
        Ptm2(i) = Ptm(k);
        Pnm2(i) = Pnm(k);
    elseif (k >= 49) && (k <= 96)
        i = k+mod(k,2);
        Ptm2(i) = Ptm(k);
        Ptm(k) = 0;
        Pnm2(i) = Pnm(k);
        Pnm(k) = 0;
    elseif (k >= 97) && (k <= 232)
        i = k + 3 - mod((k-1),4);
        Ptm2(i) = Ptm(k);
        Ptm(k) = 0;
        Pnm2(i) = Pnm(k);
        Pnm(k) = 0;
    end
end

%Niet enkel voor plots! ook voor global threshold!
i = 1;
n = 1;
a= [];
b= [];
c= [];
d= [];
for k = 1:256
   if Ptm2(k) ~= 0
       a(i) = Ptm2(k);
       b(i) = k;
       i = i+1;
   end
   if Pnm2(k) ~= 0
       c(n) = Pnm2(k);
       d(n) = k;
       n = n+1;
   end
end
figure();
hold on
plot(P);
plot(b,a, 'o');
plot(d,c, 'x');
semilogx(tq);
hold off

%Individual masking thresholds
SF = [];
SFn = [];
Ttm = [];
Tnm = [];
for i = 1:N/2
    for n = 1:N/2
        delt = FreqToBark(i*86) - FreqToBark(n*86);
        if (delt >= -3) && (delt < -1)
            SF(i,n) = 17*delt - 0.4*Ptm2(n) + 11;
            SFn(i,n) = 17*delt - 0.4*Pnm2(n) + 11;
        elseif (delt >= -1) && (delt < 0)
            SF(i,n) = (0.4*Ptm2(n)+6)*delt;
            SFn(i,n) = (0.4*Pnm2(n)+6)*delt;
        elseif (delt >= 0) && (delt < 1)
            SF(i,n) = -17*delt;
            SFn(i,n) = -17*delt;
        elseif (delt >= 1) && (delt < 8)
            SF(i,n) = (0.15*Ptm2(n) - 17)*delt - 0.15*Ptm2(n);
            SFn(i,n) = (0.15*Pnm2(n) - 17)*delt - 0.15*Pnm2(n);
        else
            SF(i,n) = -150;
            SFn(i,n) = -150;
        end
        Ttm(i,n) = Ptm2(n) - 0.275*FreqToBark(n*86) + SF(i,n) - 6.025;
        Tnm(i,n) = Pnm2(n) - 0.175*FreqToBark(n*86) + SF(i,n) - 2.025;
    end
end
%size(Ttm)
%Ttm(50,:)

figure();
hold on;
%plot(Ttm);

plot(Ttm(:,50));
plot(Ttm(:,20));
plot(Ttm(:,100));
plot(Ttm(:,150));
plot(Ttm(:,200));

plot(b,a, 'o');
%plot(Ttm(200,:));
hold off;

figure();
hold on;
%plot(Ttm);
plot(Tnm(:,50));
plot(Tnm(:,10));
plot(Tnm(:,20));
plot(Tnm(:,30));
plot(Tnm(:,100));
plot(Tnm(:,150));
plot(d,c, 'x');
%plot(Ttm(200,:));
hold off;


%Global masking thresholds

Tg = [];
for i = 1:N/2
    SumTtm = 0;
    for l = 1:length(Ttm)
       SumTtm = SumTtm + 10^(0.1*Ttm(i,l)); 
    end
    SumTnm = 0;
    for m = 1:length(Tnm)
       SumTnm = SumTnm + 10^(0.1*Tnm(i,m)); 
    end
    Tg(i) = 10*log10(10^(0.1*Tq(i))+SumTtm+SumTnm);
end

threshold = Tg;

figure();
hold on
plot(Tg);
semilogx(tq);
hold off
end

