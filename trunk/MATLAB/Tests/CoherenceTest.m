function [ output_args ] = CoherenceTest( input_args )
%COHERENCETEST 

DCTm=dctmtx(1000);

cohm=zeros(1,10);
t=1;
for i=0:50:450
    mtx=DCTm(1:end-i,:);
    [coherence, coherence_un, coherence_sq] = Coherence(mtx);
    cohm(1,t)=coherence
    t=t+1;
end

