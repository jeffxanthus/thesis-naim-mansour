function [sparsity] = estimateSparsity(signal)
%ESTIMATESPARSITY 
% Author: Naim Mansour
%This is an approximation! -TO BE EXAMINED

f=dct(signal);
maxF=max(f);
sparsity=length(find(f>0.05*maxF));
plot(f)
end

