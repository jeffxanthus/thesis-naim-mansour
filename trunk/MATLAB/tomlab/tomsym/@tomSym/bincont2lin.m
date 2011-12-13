function constraints = bincont2lin(prod, bin, cont, x_L, x_U)
% tomSym/bincont2lin - Convert product to linear constraints.
%
% constraints = bincont2lin(prod, bin, cont, x_U, x_L) creates constraints
% when modeling with binary variables which are multiplied by integer or 
% continuous variables. This is an efficient way to get rid off quadratic
% objectives or constraints. 
%
% These inequalities created:
%
% prod <= cont
% prod >= cont - x_U .* (1 - bin)
% prod <= x_U .* bin
% prod >= x_L .* bin
% cont >= x_L
% cont <= x_U
% bin  >= 0
% bin  <= 1
%
% This means that prod will always equal bin * cont
%
% INPUTS
%
% prod     product variables
% bin      binary (0-1) variables
% cont     continuous/integer variables
% x_L      lower bound on continuous variable
% x_U      upper bound on continuous variable
%
% OUTPUT
%
% constraints  A set of constraints to be included in the problem
%              definition.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

if isnumeric(x_L) && isnumeric(x_U) && any(x_L(:)>x_U(:))
    error('Upper limit must be higher than lower limit');
end

if ~isint(bin)
    error('Binary variable must be created as integer.');
end

constraints = {
    prod <= cont
    prod >= cont - x_U .* (1 - bin)
    prod <= x_U .* bin
    prod >= x_L .* bin
    cont >= x_L
    cont <= x_U
    bin  >= 0
    bin  <= 1
    };
