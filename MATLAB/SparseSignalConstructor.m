function [x,f] = SparseSignalConstructor(sparsity,length)
%SPARSESIGNALCONSTRUCTOR 

x=zeros(length,1);
for i=1:sparsity
    t=(length)*(exp(-rand()*4));
    x(round(t),1)=5;
end
f=dct(x);
subplot(2,1,1);plot(x);
subplot(2,1,2);plot(f);
end

                                                                                                                                                                                              