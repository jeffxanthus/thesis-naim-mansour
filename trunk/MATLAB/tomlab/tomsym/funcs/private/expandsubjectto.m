function fun = expandsubjectto(fun,ix,vv)
% expandsubjectto - moves out subjectTo calls, and expands symbols
%
% fun = expandsubjectto(fun,x,v) re-writes any subject-to calls inside fun
% to use one unknown per element in v. This is needed for subsVec and fsum
% to work as expected.
% It is assumed that v = 1:N.

[fv,cv,x0v] = extractConstraints(fun);
fn = fieldnames(x0v);
if ~isempty(fn)
    % Contained subjectTo - expand unknowns
    xsb = struct;
    for k=1:length(fn)
        [mm,nn] = size(x0v.(fn{k}));
        xx = tom([],length(vv),mm*nn);
        xsb.(fn{k}) = reshape(submatrix(xx,ix,':'),mm,nn);
        x0.(char(xx)) = repmat(vec(x0v.(fn{k}))',length(vv),1);
    end
    fv = subststruct(fv, xsb);
    for k=1:length(cv)
        cv{k} = subststruct(cv{k}, xsb);
    end
    fun = subjectTo(fv, x0, cv{:});
end
