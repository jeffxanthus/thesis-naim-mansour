function [f,c] = rewriteV(f,c)
% rewriteV - Rewrite an optimization problem to avoid sharp corners.
%
% [f,c] = rewriteV(f,c) rewrites an optimization problem, characterized by
% an objective function f, and a constraint set c, so as to avoid sharp
% corners in the objective funcition and constraints.
%
% f must be a scalar tomSym object.
% c must ba a cell array of tomSym constraints.
%
% Optimization alogrithms typically search for points where the derivative
% of f is zero. If the objective function contains functions like "abs"
% that have sharp corners, then such algorithms may fail to converge, or
% converge very slowly. RewriteV attemtps to rewrite the problem introducing
% extra unknowns and constraints, and eliminating the sharp corners.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-06 by rutquist for TOMLAB release 7.7

if nargin<2 || isempty(c)
    c = {};
end

if ~iscell(c);
    c = {c};
end

c = flattencellarray(c);

if isa(f,'tomSym')
    [f,cn2] = rV(f,length(f.s));
    while(~isempty(cn2))
        c = [c, cn2]; %#ok
        [f, cn2] = rV(f,length(f.s));
    end
    c = [c, cn2]; %#ok
end

i=1;
while i <= length(c)
    if tomCmp(c{i},'le')
        if c{i}.s(end).a(1) > 0 && strcmp(c{i}.s(c{i}.s(end).a(1)).op,'max')
            o1 = operand(1,c{i});
            o2 = operand(2,c{i});
            if length(o1.s(end).a)==3 && o1.s(end).a(2) < 0 && isempty(o1.d{-o1.s(end).a(2)})
                % max(a,[],dim)
                a = operand(1,o1);
                if operand(3,o1)==1
                    c{i} = a<=repmat(o2,size(a,1),1);
                else
                    c{i} = a<=repmat(o2,1,size(a,2));
                end
            else
                % max(a,b)
                c = [c(1:i-1), {operand(1,o1)<=o2, operand(2,o1)<=o2}, c(i+1:end)];
            end
            continue
        end
        if c{i}.s(end).a(1) > 0 && strcmp(c{i}.s(c{i}.s(end).a(1)).op,'abs')
            a = operand(1,operand(1,c{i}));
            o2 = operand(2,c{i});
            c = [c(1:i-1) { -o2 <= a, a <= o2 } c(i+1:end)];
            continue
        end
        [c{i}, cn2] = rV(c{i},c{i}.s(end).a(1));
        while(~isempty(cn2))
            c = [c, cn2]; %#ok
            [c{i}, cn2] = rV(c{i},c{i}.s(end).a(1));
        end
        c = [c, cn2]; %#ok
    elseif tomCmp(c{i},'ge')
        if c{i}.s(end).a(2) > 0 &&  strcmp(c{i}.s(c{i}.s(end).a(2)).op,'max')
            o2 = operand(1,c{i});
            o1 = operand(2,c{i});
            if length(o1.s(end).a)==3 && o1.s(end).a(2) < 0 && isempty(o1.d{-o1.s(end).a(2)})
                % max(a,[],dim)
                a = operand(1,o1);
                if operand(3,o1)==1
                    c{i} = a<=repmat(o2,size(a,1),1);
                else
                    c{i} = a<=repmat(o2,1,size(a,2));
                end
            else
                % max(a,b)
                c = [c(1:i-1), {operand(1,o1)<=o2, operand(2,o1)<=o2}, c(i+1:end)];
            end
            continue
        end
        if c{i}.s(end).a(2) > 0 &&  strcmp(c{i}.s(c{i}.s(end).a(2)).op,'abs')
            a = operand(1,operand(2,c{i}));
            o2 = operand(1,c{i});
            c = [c(1:i-1) { -o2 <= a, a <= o2 } c(i+1:end)];
            continue
        end
        [c{i}, cn2] = rV(c{i},c{i}.s(end).a(2));
        while(~isempty(cn2))
            c = [c, cn2]; %#ok
            [c{i}, cn2] = rV(c{i},c{i}.s(end).a(2));
        end
        c = [c, cn2]; %#ok
    end
    i = i+1;
end

function [f,cn] = rV(f,ix)
cn = {};
if ix<0
    return
end
dr = zeros(1,ix);
dr(ix) = 1;
for i=ix:-1:1
    if dr(i)==0 || isnan(dr(i))
        for k=1:length(f.s(i).a)
            if f.s(i).a(k) > 0
                dr(f.s(i).a(k)) = NaN;
            end
        end
        continue;
    end
    switch(f.s(i).op)
        case 'abs'
            if dr(i) > 0 && f.s(i).a(1) > 0
                z = tom([],f.s(i).sz1,f.s(i).sz2);
                x1 = subsymb(f.s(i).a(1),f);
                f = substsubsymb(f,i,z);
                cn = {z>=x1, z>=-x1};
                return
            end
        case 'max'
            if dr(i) > 0
                if length(f.s(i).a)==3 && f.s(i).a(2) < 0 && isempty(f.d{-f.s(i).a(2)})
                    if f.s(i).a(1) > 0 && strcmp(f.s(f.s(i).a(1)).op,'abs')
                        % max(abs(a),[],dim)
                        z = tom([],f.s(i).sz1,f.s(i).sz2);
                        x1 = subsymb(f.s(f.s(i).a(1)).a(1),f);
                        dim = f.d{-f.s(i).a(3)};
                        f = substsubsymb(f,i,z);
                        if dim==1
                            cn = {repmat(z,size(x1,1),1)>=x1,-repmat(z,size(x1,1),1)<=x1};
                        else
                            cn = {repmat(z,1,size(x1,2))>=x1,-repmat(z,1,size(x1,2))<=x1};
                        end
                        return
                    else
                        % max(a,[],dim)
                        z = tom([],f.s(i).sz1,f.s(i).sz2);
                        x1 = subsymb(f.s(i).a(1),f);
                        dim = f.d{-f.s(i).a(3)};
                        f = substsubsymb(f,i,z);
                        if dim==1
                            cn = {repmat(z,size(x1,1),1)>=x1};
                        else
                            cn = {repmat(z,1,size(x1,2))>=x1};
                        end
                        return
                   end
                elseif length(f.s(i).a)==2
                    % max(a,b)
                    z = tom([],f.s(i).sz1,f.s(i).sz2);
                    x1 = subsymb(f.s(i).a(1),f);
                    x2 = subsymb(f.s(i).a(2),f);
                    f = substsubsymb(f,i,z);
                    cn = {z>=x1, z>=x2};
                    return
                else
                    error('tomSym internal error: wrong number of args to max.');
                end
            end
        case 'uminus'
            dr = drset(dr, f.s(i).a(1), -dr(i));
        case {'uplus','sum','vec','repmat','submatrix','lookup'}
            dr = drset(dr, f.s(i).a(1), dr(i));
        case 'minus'
            dr = drset(dr, f.s(i).a(1), dr(i));
            dr = drset(dr, f.s(i).a(2), -dr(i));
        case {'plus','smplus'}
            dr = drset(dr, f.s(i).a(1), dr(i));
            dr = drset(dr, f.s(i).a(2), dr(i));
        case {'horzcat', 'vertcat'}
            for k=1:length(f.s(i).a)
                dr = drset(dr, f.s(i).a(k), dr(i));
            end
        case {'smtimes','times','mtimes'}
            if f.s(i).a(1) < 0 && ~any(diff(vec(f.d{-f.s(i).a(1)})))
                dr = drset(dr, f.s(i).a(2), dr(i)*f.d{-f.s(i).a(1)}(1));
            elseif f.s(i).a(2) < 0 && ~any(diff(vec(f.d{-f.s(i).a(2)})))
                dr = drset(dr, f.s(i).a(1), dr(i)*f.d{-f.s(i).a(2)}(1));
            end
    end
end

function dr = drset(dr, i, val)
if i > 0
    if dr(i)*val >= 0
        dr(i) = val;
    else
        dr(i) = NaN;
    end
end


