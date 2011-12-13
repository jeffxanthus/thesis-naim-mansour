function [f,c,x0] = extractConstraints(f,c,x0)
% extrectContraitns - Move constraints from subjectTo into constraint set.
%
% [f,c] = extractConstraints(f,c) modifies f and c to not contain any 
% calls to subjectTo, and appends all constraints found in subjectTo to c.
%
% This function should typically be run before converting a problem into
% mcode or attempting to compute any derivatives.
%
% See also: subjectTo

if nargin<2
    c = {};
end

if nargin<3
    x0 = struct;
end

if isa(f,'tomCmplx')
    [fr,c,x0] = extractConstraints(tomSym(real(f)),c,x0);
    [fi,c,x0] = extractConstraints(tomSym(imag(f)),c,x0);
    f = tomCmplx(fr,fi);
else
    [f,c,x0] = extractConstraints(tomSym(f),c,x0);  % Call overloaded function
end