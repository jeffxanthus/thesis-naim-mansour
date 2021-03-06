function x = tom(label, m, n, varargin)
% tom - Generate a tomSym symbol.
%
% -  x = tom creates a scalar tomSym symbol with an automatic name.
% -  x = tom(label) creates a scalar symbol with the provided name.
% -  x = tom(label,m,n) creates a m-by-n matrix symbol.
% -  x = tom([],m,n) creates a matrix symbol with an automatic name.
% -  x = tom(label,m,n,'int') creates an integer matrix symbol.
% -  x = tom(label,m,n,'symmetric') creates a symmetric matrix symbol.
%
% Because constructs like "x = tom('x')" are very common, there is the
% shorthand notation "toms x".
%
% See also: toms, tomSym

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-14 by rutquist for TOMLAB release 7.7

if ~(isempty(label) || ischar(label))
    if nargin>=3
        varargin = {n, varargin{:}};
    end
    if nargin>=2
        n = m;
    end
  m = label;
  label = [];
  warning('tomSym:MissingLabel','Missing label when creating a new tomSym.');
end

persistent idcount
if isempty(idcount)
    idcount = 1;
end
if nargin == 0 || isempty(label)
    label = ['symb' num2str(idcount)];
    idcount = idcount+1;
end

if(nargin<2)
    m = 1;
end
if(nargin<3)
    if numel(m)==2
        % Example: tom('x',size(x0));
        n=m(2);
        m=m(1);
    else
        n = 1;
    end
end

if ~isvarname(label)
    error(['Illegal tomSym symbol name ' label]);
end
if strmatch('tempC',{label})
    error('The symbol name "tempC" is reserved for autogenerated expressions.');
end
if strcmp('tempD',label)
    error('The symbol name "tempD" is reserved for autogenerated data.');
end

intflag = false;
symflag = false;

for i=1:length(varargin)
    vi = varargin{i};
    if isempty(vi)
        continue
    end
    if ~ischar(vi)
        error('Illegal flag');
    end
    if vi(1) == '-'
        vi = vi(2:end);
    end
    if strcmp(vi,'int') || strcmp(vi,'integer')
        intflag = 'integer';
    elseif strcmp(vi,'symmetric')
        symflag = true;
    else
        error(['Unknown flag: ' vi]);
    end
end

if symflag
    if m~=n
        error('Symmetric matrix must be square');
    end
    x = setSymmetric(tomSym(mfilename, n*(n+1)/2, 1, label, intflag));
else
    x = tomSym(mfilename, m, n, label, intflag);
end
