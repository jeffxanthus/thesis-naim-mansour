function y=ifThenElse(cmp1,op,cmp2,yTrue,yFalse,s)
% ifThenElse - Smoothened if/then/else.
%
% y = ifThenElse(cmp1, op, cmp2, yTrue, yFalse, s) replicates the behavior
% of the C-language construction y = ( condition ? yTrue : yFalse ) but
% with a smoothing sigmoid function, scaled by s. (Setting s=0 gives
% discontinuous "standard" if-then-else.)
%
% The operand string op must be the name of one of the comparison operators
% ('eq', 'lt', 'le', 'gt' or 'ge')
%
% For tomSym conditions, the shorened syntax is:
% y = ifThenElse(cond, yTrue, yFalse, s)
%
% The scaling factor s is optional. The default is 0.
%
% Example: The absolute value of x can be written as 
%          ifThenElse(x,'gt',0,x,-x)
%          ifThenElse(x>0,x,-x) % only works if x is a tomSym

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if islogical(cmp1)
    error('ifThenElse called using tomSym syntax, but without any tomSym object.');
end

if nargin<6 || all(s(:)<=0)
    if numel(op)==1
        if feval(op,cmp1,cmp2)
            y = yTrue;
        else
            y = yFalse;
        end
    else
        tf = double(feval(op,cmp1,cmp2));
        y = tf.*yTrue + (1-tf).*yFalse;
    end
else
    switch op
        case {'le','lt'}
            t = cmp2-cmp1;
        case {'gt','ge'}
            t = cmp1-cmp2;
        case {'eq','ne'}
            error('Not implemented yet.') % TODO
        otherwise
            error('Illegal operator for ifThenElse');
    end
    tf = 1./(1+exp(-t./s));
    y = tf.*yTrue + (1-tf).*yFalse;
end
