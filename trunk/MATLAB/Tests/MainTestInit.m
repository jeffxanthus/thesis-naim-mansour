function [SNRM1 SNRORIG1 SNRM2 SNRORIG2 SNRM3 SNRORIG3] = MainTestInit()
%MAINTESTINIT


% [SNRM2 SNRORIG2]=OMPBPTest(1,'BeethoP5.wav');
% [SNRM3 SNRORIG3]=OMPBPTest(1,'Forsake2.wav');
% [SNRM2 SNRORIG2]=OMPBPTest(1,'BachHymn.wav');
[SNRM1 SNRORIG1]=OMPBPTest(1,'bach_partita.wav');
[SNRM3 SNRORIG3]=OMPBPTest(1,'Secret_1.wav');

% [RM4 SNRM4]=AxBeTest(3,'Secret_2.wav');

end


