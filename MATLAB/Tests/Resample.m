function [] = Resample(fileName,fsNew)
%RESAMPLE 

[data,fs,noBits]=wavread(fileName);
data=resample(data,fsNew,fs);
wavwrite(data,fsNew,noBits,[fileName '2']);
end

