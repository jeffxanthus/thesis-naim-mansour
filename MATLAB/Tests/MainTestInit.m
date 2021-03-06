function [SNRM1 SNRORIG1 ODG1 ODGORIG1 Rsig1 Rsam1 SNRM3 SNRM2 SNRORIG2 ODG2 ODGORIG2 Rsig2 Rsam2] = MainTestInit()
%MAINTESTINIT


% [SNRM2 SNRORIG2]=OMPBPTest(1,'BeethoP5.wav');
% [SNRM3 SNRORIG3]=OMPBPTest(1,'Forsake2.wav');
% [SNRM2 SNRORIG2]=OMPBPTest(1,'BachHymn.wav');

[SNRM1 SNRORIG1 ODG1 ODGORIG1 Rsig1 Rsam1]=OMPBPTest(2,'bach_partita.wav',30,[0.76  0.62  0.55],[0.86 0.77 0.68],[3]);
save('trueL0unboundedPatchedBachPartitaAxBe.mat', 'SNRM1', 'SNRORIG1', 'ODG1', 'ODGORIG1','Rsig1', 'Rsam1');
[SNRM2 SNRORIG2 ODG2 ODGORIG2 Rsig2 Rsam2]=OMPBPTest(2,'Secret_1.wav',30,[0.38  0.31  0.26],[0.37 0.31 0.25],[3]);
% [SNRM2 SNRORIG2 ODG2 ODGORIG2 Rsig2 Rsam2]=OMPBPTest(2,'FixYou.wav',30,[0.47 0.40 0.32],[0.43 0.35 0.28],[1 2]);
save('trueL1unboundedPatchedSecret1AxBe.mat', 'SNRM2', 'SNRORIG2', 'ODG2', 'ODGORIG2','Rsig2','Rsam2');
% [SNRM3 SNRORIG3 ODG3 ODGORIG3 Rsig3 Rsam3]=OMPBPTest(2,'FixYou.wav',30,[0.47 0.40 0.32],[0.43 0.35 0.28],[3]);
% save('trueL1unboundedPatchedFixYouAxBe.mat', 'SNRM3', 'SNRORIG3', 'ODG3', 'ODGORIG3','Rsig3','Rsam3');
% [SNRM4 SNRORIG4 ODG4 ODGORIG4 Rsig4 Rorig4]=OMPBPTest(1,'bach_partita.wav',30,[0.76  0.62  0.55],[0.86 0.77 0.68],[4]);
% save('IRL1unboundedPatchedBachPartitaAxBe.mat', 'SNRM4', 'SNRORIG4', 'ODG4', 'ODGORIG4');

% [a5 b5 SNRM5 SNRORIG5]=OMPBPTest(1,'Secret_2.wav',30);
% [a6 b6 SNRM6 SNRORIG6]=OMPBPTest(1,'FolkMusic.wav',30);
% [a7 b7 SNRM7 SNRORIG7]=OMPBPTest(1,'FolkMusic.wav',64);

% [SNRM3 SNRORIG3]=OMPBPTest(1,'Secret_1.wav');

% [RM4 SNRM4]=AxBeTest(3,'Secret_2.wav');

end


