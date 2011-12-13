function [f,c,x0] = extractConstraints(f,c,x0)
% extrectContraitns - Move constraints from subjectTo into constraint set.
%
% [f,c,x0] = extractConstraints(f,c,x0) modifies f and c to not contain any 
% calls to subjectTo, and appends all constraints found in subjectTo to c,
% and their initial guesses to the struct x0.
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

if isa(f,'tomSym');
    % Extract constraints from f
    
    done = false;
    while ~done
        done = true;
        for i=1:length(f.s)
            if strcmp(f.s(i).op,'subjectTo')
                for k=3:length(f.s(i).a)
                    c{end+1} = subsymb(f.s(i).a(k),f);
                end
                x0k = f.d{-f.s(i).a(2)};
                fn = fieldnames(x0k);
                for k=1:length(fn);
                    x0.(fn{k}) = x0k.(fn{k});
                end
                f = substsubsymb(f,i,subsymb(f.s(i).a(1),f));
                done = false;
                break;
            end
        end
    end
end

if ~isempty(c)
    % Recursively extract constraints from c.
    if ~iscell(c);
        c = { c };
    end
    c = c(:)';
    
    i = 1;
    while i<=length(c)
        if iscell(c{i})
            c = {c{1:i-1}, c{i}{:} c{i+1:end}};
            continue
        end
        if islogical(c{i})
            i = i+1;
            continue
        end
        [ci,cr,x0] = extractConstraints(c{i},{},x0);
        c{i} = ci;
        for k=1:length(cr)
            % Check that the same constraint is not added more than once.
            iseq = false;
            for m = 1:length(c);
                if isequal(c{m},cr{k})
                    iseq = true;
                    break;
                end
            end
            if ~iseq
                c{end+1} = cr{k};
            end
        end
        i = i+1;
    end
end
