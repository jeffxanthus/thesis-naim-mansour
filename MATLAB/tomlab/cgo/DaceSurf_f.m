% function f = DaceSurf_f(x, Prob)
%
% DACE Surface Evaluation

function f = DaceSurf_f(x, Prob)

f = cgolib(205, Prob.EGO.CGOLIB.daceid(Prob.EGO.CGOLIB.TRANSFORM+1), x);

% MODIFICATION LOG
%
% 080112 hkh  Written
