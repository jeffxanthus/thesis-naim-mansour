function [A, B, idx] = sym2eig(lambda, c)
% sym2eig Convert symbolic constraints into an eigenvalue problem
%
% [A, B, S, idx] = sym2eig(lambda, c) converts the symbolic constraints c
% into a generalized eigenvalue problem A*x = lambda*B*x
%
% Inptus:
% lambda must be a tomSym symbol.
% c must be a cell array of tomSym equations compatible with the
% eigenvalue problem formulation.
%
% Outputs:
% A and B will be numeric matrices.
% idx will be a struct containing the indexes of the symbols (other than
% lambda). 
%
% If c contains any equations of the type x(i) == 0, then the corresponding
% elements of idx will be 0, there will be no corresponding column in the
% matrices A and B.

% Expand list of constraints. (Code copied from sym2prob and slightly
% modified).

if ~iscell(c)
    c = {c};
end

i = 1;
while i<=length(c)
    if iscell(c{i})
        c = {c{1:i-1}, c{i}{:} c{i+1:end}};
        continue
    end
    if tomCmp(c{i},'vertcat') || tomCmp(c{i},'horzcat') || ...
            tomCmp(c{i},'cellarray')
        a = operands(c{i});
        c = {c{1:i-1}, a{:}, c{i+1:end}};
        continue
    end
    if isa(c{i}, 'tomArray')
        c{i} = vec(c{i});
        if tomCmp(c{i},'vec')
            c{i} = operand(1,c{i});
        end
    end
    if isempty(c{i})
        c = {c{1:i-1} c{i+1:end}};
        continue
    end
    if ~isa(c{i}, 'tomSym')
        if ~isa(c{i},'logical')
            error(['Illegal constraint type: "' class(c{i}) ...
                '" - Constraints must be either tomSym or logical.']);
        end
        if ~all(all(c{i}))
            error('Constant constraint evaluated to "false".');
        end
        warning('tomSym:ConstConstraint','Constant constraint removed.');
        c = {c{1:i-1}, c{i+1:end}};
        continue
    end
    if strcmp(operator(c{i}),'positiveSemidefinite')
        error('Constraint incompatible with sym2eig');
    end
    if tomCmp(c{i},'dcollocate')
        error('Missing call to docollocate?');
    end
    if tomCmp(c{i},'ctranspose')
        c{i} = operand(1,c{i});
    end
    if tomCmp(c{i},'complementary')
        error('Constraint incompatible with sym2eig');
    end
    if tomCmp(c{i},'and')
        % "and" has the same effect as listing the operands as separate
        % constraints
        c = {c{1:i-1}, operand(1,c{i}), operand(2,c{i}), c{i+1:end}};
        continue
    end
    if isempty(strmatch(operator(c{i}),{'eq'}))
        error('Only equality constraints are allowed.');
    end
    if isa(operand(1,c{i}),'tomSym') && ...
            ~isempty(strmatch(operator(operand(1,c{i})),{'eq','le','ge'}))
        % Expand (a == b) == c
        c = {c{1:i-1}, operand(1,c{i}), ...
            feval(operator(c{i}), operand(2,operand(1,c{i})), operand(2,c{i})), ...
            c{i+1:end}};
    elseif isa(operand(2,c{i}),'tomSym') && ...
            ~isempty(strmatch(operator(operand(2,c{i})),{'eq','le','ge'}))
        % Expand a == b == c
        c = {c{1:i-1}, ...
            feval(operator(c{i}), operand(1,c{i}), operand(1,operand(2,c{i}))), ...
            operand(2,c{i}), c{i+1:end}};
    else
        i = i+1; % Proceed to next constraint.
    end
end

% Compile list of symbols, excluding params.
symbs = symbols(tomSym({lambda,c}),'struct');

symbs = rmfield(symbs,char(lambda));

nx = 0;
idx = struct;
symbs = orderfields(symbs);
sl = fieldnames(symbs);
x = cell(size(sl));
for i=1:length(sl)
    x{i} = symbs.(sl{i});
    idx.(sl{i}) = reshape(nx+1:nx+numel(x{i}),size(x{i}));
    nx = nx + numel(x{i});
    if operand(2,x{i})
        error('Integer variables not allowed with sym2eig');
    end
    x{i} = vec(x{i});
end
x = vertcat(x{:});

A = zeros(length(x),length(x));
B = zeros(length(x),length(x));
isz = false(1,length(x));
n = 0;
for i=1:length(c)
    % c should be only equality constraints by now.
    ci = operand(2,c{i}) - operand(1,c{i});
    ni = numel(c{i});
    Ai = derivative(subs(ci,lambda==0),x);
    Bi = -derivative(derivative(ci,lambda),x);
    if ~( all(vec(constpart(ci))==0) && isnumeric(Ai) && isnumeric(Bi) )
        error(['Constraint not linear: ' mcodestr(c{i})])
    end
    if nnz(Bi)==0 && all(sum(Ai~=0,2)==1)
        isz = sum([isz;Ai~=0])~=0;      
    else
        A(n+1:n+ni,:) = Ai;
        B(n+1:n+ni,:) = Bi;
        n = n+ni;
    end
end

if any(isz)
    newix = cumsum(isz==0);
    newix(isz==1) = 0;
    sl = fieldnames(idx);
    for i=1:length(sl)
        idx.(sl{i}) = newix(idx.(sl{i}));
    end
    A = A(1:n,isz==0);
    B = B(1:n,isz==0);
end

if n~=length(isz)-nnz(isz)
    warning('sym2eig:nonsquare', ...
        ['Expected ' num2str(length(x)) ' constraints, found ' num2str(n) '.']);
    A = A(1:n,:);
    B = B(1:n,:);
end

[scratch,ix] = sortrows(B~=0); %#ok
ix = flipud(ix);
A = A(ix,:);
B = B(ix,:);

if B(1,1)<0
    % Prefer positive values in B. (B should be pos semidefinite.)
    A = -A;
    B = -B;
end

if ~all(vec(B==B'))
    warnign('sym2eig:nonsymB','B is not symmetric.');
end
