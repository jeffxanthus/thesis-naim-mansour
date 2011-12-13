function Prob=setParameter(Prob,varargin)
% setParameter - Give a fixed value to a tomSym symbol in a Prob
%
% y = Prob=setParameter(Prob,symbol1,value1,...)
% y = Prob=setParameter(Prob,'symbol1',value1,...)
% y = Prob=setParameter(Prob,struct('symbol1,value1,...))
% y = Prob=setParameter(Prob,symbol1==value1,...) 
%
% setParameter is used to give a fixed value to one or more of the tomSym
% symbols used in a problem. This is useful for quickly re-solving the same
% problem, using different values of a parameter without generating new
% m-code each time.
%
% Re-solving a problem can be an effective way of obtaining convergence for
% strongly non-linear problems. - Using a sequence of problems such that
% the first one is nearly linear, re-solve the problem using the previous
% solution as starting guess for the next one.
%
% Example:
% toms a
% alist = [...]; % List of a-values to solve for
% ...
% Prob = sym2prob(...);
% for i=1:length(alist)
%    Prob = setParameter(Prob,a==alist(i));
%    Result = tomRun('solver',Prob);
%    Prob.x_0 = Result.x_k; % Use solution as starting guess for next round
% end
% solution = getSolution(Result);
%
% See also: sym2prob tomRun subs

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if length(varargin)==1 && isstruct(varargin{1})
    s = varargin{1};
elseif iscell(varargin{1}) || tomCmp(varargin{1},'eq')
    s = tom2struct(varargin);
else
    s = struct;
    for i=1:2:length(varargin)-1;
        symb = varargin{i};
        if isa(symb, 'tomSym')
            symb = operand(1,symb);
        end
        if ~ischar(symb)
            error(['Expected a symbol name. Got ' class(symb) '.']);
        end
        s.(symb) = varargin{i+1};
    end
end

fn = fieldnames(s);
for i=1:length(fn)
    Prob.x_L(vec(Prob.tomSym.idx.(fn{i}))) = vec(s.(fn{i}));
    Prob.x_U(vec(Prob.tomSym.idx.(fn{i}))) = vec(s.(fn{i}));
    Prob.x_0(vec(Prob.tomSym.idx.(fn{i}))) = vec(s.(fn{i}));
end
