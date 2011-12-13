function Prob = autoScaleProb(Prob)
% autoScaleProb - Automatically set scaling for a scaleable tomSym problem
%
% Prob = autoScaleProb(Prob) changes the scaling used in Prob to values
% determined auomatically by analyzing the problem limits and derivatives.
%
% See also: scaleProb

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2011 by Tomlab Optimization Inc.
% Last modified 2011-02-28 by rutquist for TOMLAB release 7.7

if ~isfield(Prob,'tomSym')
    error('The Prob struct does not seem to originate from a tomSym problem.');
end

if ~isfield(Prob.tomSym,'xScale')
    error('The problem does not seem to use scaling. (Try setting options.scale)');
end

% ------- Find proper xScale and guessing range ------
gL = Prob.x_L;
gU = Prob.x_U;
xs = gU-gL;
x0 = Prob.x_0;
    
idx = Prob.tomSym.idx;
fn = fieldnames(idx);
for i=1:length(fn)
    ix = vec(idx.(fn{i}))';

    if any(isinf(xs(ix)))
        if any(~isinf(xs(ix)))
            % If limits exist for any member of the set,
            % then use that for scaling the entire set.
            xs(ix(isinf(xs(ix)))) = max(xs(ix(~isinf(xs(ix)))));
        else
            % No member of the set has limits.
            % Use starting guess for scaling.
            xs(ix) = max(abs(x0(ix)));
        end
    end

    if any(xs(ix)>0)
        % Don't let the range become too extreme - stay within one order of
        % magnitude of the mean.
        mxs = mean(xs(ix(xs(ix)>0)));
        xs(ix) = max(0.1*mxs,min(10*mxs,xs(ix)));
    else
        % Scale can't be zero - reset to 1.
        xs(ix) = 1;
    end
end

% Integer variables cannot be scaled - reset to 1.
if isfield(Prob.MIP,'IntVars')
    xs(Prob.MIP.IntVars) = 1;
end

% Fix limits for guesses
ix = isinf(gL);
gL(ix) = x0(ix) - 0.5*xs(ix);
ix = isinf(gU);
gU(ix) = x0(ix) + 0.5*xs(ix);
if isfield(Prob.MIP,'IntVars')
    gL(Prob.MIP.IntVars) = ceil(gL(Prob.MIP.IntVars));
    gU(Prob.MIP.IntVars) = floor(gL(Prob.MIP.IntVars));
end

% ------- Create a set of guesses ------
gv = zeros(length(x0),5);
for i=1:size(gv,2);
    gv(:,i) = x0 + 0.1*sin(i*(1:length(x0))').*(gU-gL);
    gv(:,i) = min(gU,max(gL,gv(:,i)));
end

% ------- Find scale for constraints ------
if isfield(Prob.tomSym,'cScale') && ~isempty(Prob.FUNCS.dc)
    % Find scale for c.
    cv = zeros(length(Prob.c_L),5);
    for i=1:size(cv,2);
        dc = feval(Prob.FUNCS.dc,gv(:,i),Prob);
        cv(:,i) = sqrt(sum((repmat(xs',size(dc,1),1).*dc).^2,2));
    end
    cv(cv<=0)=1;
    scv = exp(mean(log(cv),2));
    %scv(scv==0)=1;
    sidx = Prob.tomSym.lambdaIdxNonlin;
    eidx = [Prob.tomSym.lambdaIdxNonlin(2:end)-1; size(cv,1)];
    cs = ones(size(scv));
    for i=1:length(sidx)
        cs(sidx(i):eidx(i)) = 1./exp(mean(log(scv(sidx(i):eidx(i)))));
    end
    % TODO: Group constraints
    cs(isnan(cs)) = 1;
else
    cs = [];
end

if isfield(Prob.tomSym,'fScale') && ~isempty(Prob.FUNCS.g) && ~isempty(cs)
    % Find scale for f. Note that actually rescaling f will change the
    % objective value, as displayed by the solver.
    fv = zeros(5,1);
    for i=1:length(fv);
        fv(i) = norm(xs.*feval(Prob.FUNCS.g,gv(:,i),Prob),2);
    end
    fv(fv<=0) = 1;
    sfv = exp(mean(log(fv)));
    fs = 1./sfv;
else
    fs = [];
end

if ~isempty(Prob.A)
    % Find scale for A.
    asc = ones(size(Prob.A,1),1);
    for i=1:length(asc)
        asc(i) = norm(xs'.*Prob.A(i,:),2);
    end
    asc(asc==0) = 1;
    sidx = Prob.tomSym.lambdaIdxLin;
    eidx = [Prob.tomSym.lambdaIdxLin(2:end)-1; size(Prob.A,1)];
    as = ones(size(asc));
    for i=1:length(sidx)
        as(sidx(i):eidx(i)) = 1./exp(mean(log(asc(sidx(i):eidx(i)))));
    end
    as(isnan(as)) = 1;
else
    as = [];
end

% Actually: Let's get rid of fs. The only thing that matters is that c, a
% and f are the same order of magnitude.
if ~isempty(fs)
    cs = cs./fs;
    as = as./fs;
    fs = [];
end

Prob = scaleProb(Prob, xs, cs, as, fs);
