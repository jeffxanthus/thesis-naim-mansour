function [SNRM1 SNRORIG1 ODG1 ODGORIG1 Rsig1 Rsam1 SNRM2 SNRORIG2 ODG2 ODGORIG2 Rsig2 Rsam2 SNRM3 SNRORIG3 ODG3 ODGORIG3 Rsig3 Rsam3] = MainTestInit()
%MAINTESTINIT


% [SNRM2 SNRORIG2]=OMPBPTest(1,'BeethoP5.wav');
% [SNRM3 SNRORIG3]=OMPBPTest(1,'Forsake2.wav');
% [SNRM2 SNRORIG2]=OMPBPTest(1,'BachHymn.wav');
[SNRM1 SNRORIG1 ODG1 ODGORIG1 Rsig1 Rsam1]=OMPBPTest(1,'bach_partita.wav',30);
save('test1.mat', 'SNRM1', 'SNRORIG1', 'ODG1', 'ODGORIG1','Rsig', 'Rsam');
[SNRM2 SNRORIG2 ODG2 ODGORIG2 Rsig2 Rsam2]=OMPBPTest(1,'BeethoP5.wav',30);
save('test2.mat', 'SNRM2', 'SNRORIG2', 'ODG2', 'ODGORIG2');
[SNRM3 SNRORIG3 ODG3 ODGORIG3 Rsig3 Rsam3]=OMPBPTest(1,'FolkMusic.wav',30);
save('test3.mat', 'SNRM3', 'SNRORIG3', 'ODG3', 'ODGORIG3');
% [SNRM4 SNRORIG4 ODG4 ODGORIG4 Rsig4 Rsam4]=OMPBPTest(1,'Secret_1.wav',30);
% save('test4.mat', 'SNRM4', 'SNRORIG4', 'ODG4', 'ODGORIG4');

% [a5 b5 SNRM5 SNRORIG5]=OMPBPTest(1,'Secret_2.wav',30);
% [a6 b6 SNRM6 SNRORIG6]=OMPBPTest(1,'FolkMusic.wav',30);
% [a7 b7 SNRM7 SNRORIG7]=OMPBPTest(1,'FolkMusic.wav',64);

% [SNRM3 SNRORIG3]=OMPBPTest(1,'Secret_1.wav');

% [RM4 SNRM4]=AxBeTest(3,'Secret_2.wav');

end


