% expGetLa picks up alpha and lambda parameters for exponential fitting from
% problem structure Prob.ExpFit
%
% function [lambda, alpha] = expGetLa(Prob);
%
% INPUT:
%  Prob     Problem structure
%
% OUTPUT:
%  lambda   p-dimensional lambda vector (intensities)
%  alpha    p-dimensional alpha vector (weights)
%           p is stored in Prob.ExpFit.p

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Feb 17, 1998.  Last modified June 21, 1999.

function [lambda, alpha, beta] = expGetLa(Prob)

if isempty(Prob)
    lambda=[];
    alpha =[];
    beta=[];
else
    lambda=Prob.ExpFit.lambda;
    alpha =Prob.ExpFit.alpha;
    beta =Prob.ExpFit.beta;
end