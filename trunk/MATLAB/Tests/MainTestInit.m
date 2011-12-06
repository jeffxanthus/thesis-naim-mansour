function [SNRM1 SNRORIG1 SNRM2 SNRORIG2 SNRM3 SNRORIG3 SNRM4 SNRORIG4 SNRM5 SNRORIG5 SNRM6 SNRORIG6 SNRM7 SNRORIG7] = MainTestInit()
%MAINTESTINIT


% [SNRM2 SNRORIG2]=OMPBPTest(1,'BeethoP5.wav');
% [SNRM3 SNRORIG3]=OMPBPTest(1,'Forsake2.wav');
% [SNRM2 SNRORIG2]=OMPBPTest(1,'BachHymn.wav');
[a1 b1 SNRM1 SNRORIG1]=OMPBPTest(1,'bach_partita.wav',32);
[a2 b2 SNRM2 SNRORIG2]=OMPBPTest(1,'bach_partita.wav',48);
[a3 b3 SNRM3 SNRORIG3]=OMPBPTest(1,'bach_partita.wav',64);
[a4 b4 SNRM4 SNRORIG4]=OMPBPTest(1,'BeethoP5.wav',48);
[a5 b5 SNRM5 SNRORIG5]=OMPBPTest(1,'BeethoP5.wav',64);
[a6 b6 SNRM6 SNRORIG6]=OMPBPTest(1,'FolkMusic.wav',48);
[a7 b7 SNRM7 SNRORIG7]=OMPBPTest(1,'FolkMusic.wav',64);

% [SNRM3 SNRORIG3]=OMPBPTest(1,'Secret_1.wav');

% [RM4 SNRM4]=AxBeTest(3,'Secret_2.wav');

end


