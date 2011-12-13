function o = subststruct(o,s,cs)
% tomSym/subststruct - Substitue tomSym symbols from a struct
%
% o = subststruct(o,s,cs) substitutes the symbols defined in struct into
% the symbolic expression o.
%
% If cs = true, then "constant" symbols are replaced by their numeric
% value.
%
% See also: subs

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-08-22 by rutquist for TOMLAB release 7.7

if ~isa(o, 'tomSym')
    return
end

if nargin<3
    cs = false;
end

slst = [];

for i=1:length(o.s)
    if strcmp(o.s(i).op,'subjectTo')
        if isempty(slst)
            slst = symbols(o,'all');
        end
        % Create new unknowns to avoid conflicts other substitutions
        a2 = o.d{-o.s(i).a(2)};
        fn = fieldnames(a2);
        for j=1:length(fn)
            s.(fn{j}) = tomreplacement(slst.(fn{j}));
        end
        % We also need to make subsitutions inside initial guess,
        % and change field names in struct.
        a2n = struct;
        for j=1:length(fn)
            if isa(a2.(fn{j}),'tomSym')
                a2n.(char(s.(fn{j}))) = subststruct(a2.(fn{j}),s,cs);
            else
                a2n.(char(s.(fn{j}))) = a2.(fn{j});
            end
        end
        o.d{-o.s(i).a(2)} = a2n;
    end
end

if ~cs && isempty(fieldnames(s))
    return
end

ns = cell(1,length(o.s));  % ns will contain the new subsymbols
a = false(1,length(o.s));  % a will be true if the subsymbol has changed

for i=1:length(o.s)
    if strcmp(o.s(i).op,'tom')
        xn = o.d{-o.s(i).a(1)};
        if isfield(s,xn)
            % Substitute symbol
            ns{i} = s.(xn);
            a(i) = true;
        else
            ns{i} = subsymb(i,o);
        end
    elseif strcmp(o.s(i).op,'constant')
        if cs
            ns{i} = o.d{-o.s(i).a(1)};
            a(i) = true;
        else
            ns{i} = subsymb(i,o);
        end     
    elseif strcmp(o.s(i).op,'srepmat') && o.s(i).a(1) < 0
        if cs
            ns{i} = repmat(o.d{-o.s(i).a(1)},o.d{-o.s(i).a(2)},o.d{-o.s(i).a(3)});
            a(i) = true;
        else
            ns{i} = subsymb(i,o);
        end     
    elseif ~any(a(o.s(i).a(o.s(i).a>0)))
        ns{i} = subsymb(i,o);
    else
        arg = cell(size(o.s(i).a));
        for k=1:length(arg)
            if o.s(i).a(k) < 0
                arg{k} = o.d{-o.s(i).a(k)};
            else
                arg{k} = ns{o.s(i).a(k)};
            end
        end
        ns{i} = builtin('feval',o.s(i).op,arg{:});
        a(i) = true;
    end
end

% Cleanup
o = ns{end};
