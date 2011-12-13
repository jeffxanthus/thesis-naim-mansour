function os=subs(o,varargin)
% Subs - Substitution for tomSym objects
%
% y = subs(x,symbol1,value1,...)
% y = subs(x,'symbol1',value1,...)
% y = subs(x,struct('symbol1,value1,...))
% y = subs(x,{symbol1==value1,...})
% y = subs(x,solution) % where solution was returned by tomRun
%
% If O is a cell array or a structure, then subs is applied recursively
% to all fields of the structure or array.
%
% SUBS('WORKSPACE',...) applies SUBS recursively to every variable
% in the calling workspace. 
%
% If many different values should be substituted it will be faster to
% create a .m file from the tomSym object using the "mfile" command, and
% then calling that function instead of subs.
%
% See also: tomSym/subs mfile

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-08-22 by rutquist for TOMLAB release 7.7

% TODO: Accept sub(symb,old,new) syntax.

if length(varargin)==1 && isstruct(varargin{1})
    s = varargin{1};
    if isfield(s,'Prob') && isstruct(s.Prob)
        if ~isempty(s.Prob.tomSym.phase)
            o = docollocate(s.Prob.tomSym.phase,o);
        end
       s = getSolution(s);
     end
elseif length(varargin)==1 && (iscell(varargin{1}) || tomCmp(varargin{1},'eq'))
    s = tom2struct(varargin{1});
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

if isa(o,'tomSym')
    os = subststruct(o,s,1);
elseif isa(o,'cell')
    os = cell(size(o));
    for i=1:numel(o)
        os{i} = subs(o{i},s);
    end
elseif isa(o,'struct')
    os = struct;
    fn = fieldnames(o);
    for i=1:length(fn)
        os.(fn{i}) = subs(o.(fn{i}),s);
    end        
elseif ischar(o) && strcmpi(o,'workspace') 
    % Set variables in caller workspace.
    wsvar = evalin('caller','who');
    os = 0;
    for i=1:length(wsvar)
        v = evalin('caller',wsvar{i});
        if isa(v,'tomSym')
            vs = subs(v,s);
            assignin('base',wsvar{i},vs);
            os = os+1;
        elseif isa(v,'cell') || isa(v,'struct')
            vs = subs(v,s);
            assignin('base',wsvar{i},vs);
            os = os+1;
        end
    end
else
    os = o;
end
