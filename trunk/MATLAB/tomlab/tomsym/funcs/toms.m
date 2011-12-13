function toms(varargin)
% toms - Create tomSym objects.
%
% Toms is a shorthand notation, possibly replacing several calls to 'tom'
%
% A symbol is created in the current workspace for each name listed. If a
% size is specified on the format "NxM" where N and M are integers, then
% all subsequent symbols will get that size. If the size specification
% ends with an exclamation point (as in "3x4!") then a symbolic array of
% concatenated scalar symbols is created instead of one matrix symbol.
%
% The flags "integer" (or "int") and "symmetric" are recognized. If a flag
% is encountered, then all subsequent symbols will get the properties of
% that flag.
%
% - integer:   The variable is constrained to integer values, resulting in
%              a mixed-integer problem (which requires a compatible solver.)
%
% - symmetric: The variable is symmetric, so that x' == x. This requires 
%              that the dimensions be square. An NxN symmetric matrix only 
%              contains N*(N+1)/2 unknowns, and the resulting symbolic 
%              object uses the setSymmetric function.
%
% Examples:
%
%   toms x y z
% is equivalent to
%   x = tom('x');
%   y = tom('z');
%   z = tom('z');
%
%   toms 2x3 Q 3x3 -integer R -symmetric S
% is equivalent to
%   Q = tom('Q', 2, 3);
%   R = tom('R', 3, 3, 'integer');
%   S = tom('S', 3, 3, 'integer', 'symmetric')
%
%   toms 3x1! v
% is equivalent to
%   v1 = tom('v1'); 
%   v2 = tom('v2'); 
%   v3 = tom('v3'); 
%   v = [v1;v2;v3];
%
% In the last example, with the exclamation point, the result is a vector
% containing scalar symbols. This works differently from the matrix symbols
% used in the previous examples. Mathematically this v is equivalent to 
% "toms 3x1 v", but the auto-generated code will be different. Expressions
% such as v(1)+sin(v(2))*v(3) will be more efficient, while expressions
% such as A*v will be less efficient.
%
% Note: While the "toms" shorthand is very convenient to use prototyping
% code, it is recommended to only use the longhand "tom" notation for 
% production code. The reason is that Matlab's internal compiler tries to
% guess whether a statement like "x(1)" is an index into a vector or a call
% to a function. Since it does not realize that the call to "toms" creates
% new variables is will make the wrong guess if that function (named "x" in 
% this example) is on the Matlab path when the script is loaded into
% memory. This will cause strange and unexpected results.
%
% See also: tom, tomSym, setSymmetric

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-10-05 by rutquist for TOMLAB release 7.7

y = toms_helper('tom', {}, varargin{:});

for i=1:length(y)
    if ~isempty(y{i})
        assignin('caller', varargin{i}, y{i});
    end
end
