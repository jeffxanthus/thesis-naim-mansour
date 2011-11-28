function [ H ] = perceptualWeightingMatrix( maskingThreshold)
%perceptualWeightingMatrix creates filtering matrix based on the given
%maskingThreshold.
%   Steven De Hertogh

K = length(maskingThreshold);
H = zeros(K,K);
A = 1./maskingThreshold;

figure();
plot(maskingThreshold)

%plot(A);

%Generate the filter coefficients
h = [];
for n = 1:K
   sum = 0;
   for k = 1:K
       sum = sum + sqrt(A(k))*cos(2*pi*k*n/K);
   end

   h(n) = sum/K;
end

%Generate the matrix
    %bottom half
s=0;
for i = 1:K
    n = 1;
    while n+s <= K 
      H(n+s,i) = h(n); 
      n = n+1;
    end
    s = s+1;
end
    %top half
s=0;
for n = 1:K
   i = 2;
   while i+s <= K
       H(n,i+s) = h(K-(i-2));
       i = i+1;
   end
   s=s+1;
end

end

