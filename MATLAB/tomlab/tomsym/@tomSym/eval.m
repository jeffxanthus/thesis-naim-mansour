function y = eval(o)
% tomSym/eval - Evaluate a tomSym object in the current workspace
%
% y = eval(o) evaluates the object o in the current workspace.
%
% This is useful if any of the symols used to create the object has
% changed.
%
% Example:
%    toms x y;
%    z = x*y;
%    x = 4;
%    y = 2;
%
% Now disp(z) still produces the result x*y, while eval(z) yields 8.
%
% See also: tomSym/subs

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

symbs = symbols(o);

s = struct;
for i=1:length(symbs)
    if evalin('caller',['exist(''' symbs{i} ''',''var'')'])
        s.(symbs{i}) = evalin('caller',symbs{i});
    end
end

y = subststruct(o,s);
