function display(p)
% tomArray/display - Command window display of a tomArray.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

disp(' ');
str = [inputname(1),' = tomArray('];
for i=1:length(p.sz)
    if ~isempty(p.ni)
	    str = [str p.ni{i} '='];
    end
	str = [str '1..' num2str(p.sz(i))];
    if i<length(p.sz)
        str = [str ', '];
    end
end
if isa(p.X, 'tomSym')
    disp([str ') with tomSym data:']);
    disp(' ');
    disp(['   ' char(p.X)]);
else
    disp([str ') with numeric data']);
end
disp(' ');
