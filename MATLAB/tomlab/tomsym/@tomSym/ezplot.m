function h=ezplot(varargin)
% tomSym/ezplot - Function plot
%
% ezplot(y,a,b) where y is a tomSym object containting a single scalar
% variable, plots y as a function of that variable from a to b.
%
% ezplot(x,y,a,b) creates a parametric plot.
%
% ezplot(x) plots a trajectory for the tomState or tomControl x using
% the last solution computed by ezsolve.
%
% If x and/or y contain more than one variable, ezplots tries to substitute
% the last solution obtained by ezsolve (which tomSym stores in the global
% variable lastTomSymSolution). If this results in a function of only one
% variable, then it is plotted.
%
% If the a and/or b arguments are symbolic, ezplot tries to substitute the
% last solution into them as well.
%
% If the a and b arguments are omitted, ezplot checks for a propt phase
% (stored in the last ezsolve solution or in the global varaiable proptActivePhase).
% If such a phase is found, the intervall length from the phase is used,
% and the interpolation points are marked in the plot.
%
% Example:
%
%  toms t
%  figure(1)
%  ezplot(sin(t),0,2*pi)        % plots a sine wave
%  figure(2)
%  ezplot(cos(t),sin(t),0,2*pi) % plots a circle

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-01-14 by rutquist for TOMLAB release 7.7

global lastTomSymSolution
global proptActivePhase

if isempty(varargin) || ~isa(varargin{1},'tomSym')
    error('EZPLOT expects first argument to be a symbolic object');
end

na = length(varargin);

if ~isempty(lastTomSymSolution) && ~isempty(lastTomSymSolution.Prob.tomSym.phase)
    phase = lastTomSymSolution.Prob.tomSym.phase;
elseif ~isempty(proptActivePhase)
    phase = proptActivePhase;
else
    phase = [];
end

if length(symbols(varargin{1})) > 1 || ~isempty(phase) || ...
        ((na==2 || na>=4) && length(symbols(varargin{2})) > 1) || ...
        (na==3 && isa(varargin{2},'tomSym')) || ...
        (na>=3 && isa(varargin{3},'tomSym')) || ...
        (na>=4 && isa(varargin{4},'tomSym'))
    % Substitution of solution is necessary.
    if ~isempty(lastTomSymSolution)
        varargin = subs(varargin,lastTomSymSolution);
    end
end

x = varargin{1};
if isnumeric(x) && ~isempty(phase)
    % A PROPT solution happened to be constant.
    x = tomSym(x);
    s = phase.t;
elseif length(symbols(x)) == 1
    s = symbols(x,'struct');
    s = struct2cell(s);
    s = s{1};
else
    error('Input X is not a function of a single symbol.');
end

if na>=2 && isa(varargin{2},'tomSym')
    % Y argument provided
    y = varargin{2};
    if length(symbols([vec(x); vec(y)])) ~= 1
        error('Inputs X and Y are not functions of a single symbol.');
    end
    ny = 1;
else
    ny = 0;
end

if na>=ny+3 && ~ischar(varargin{ny+2})
    % a and b were provided
    a = varargin{ny+2};
    b = varargin{ny+3};
    nab = 2;
    coll = false;
    npts = 100;
else
    if ~isempty(phase)
        a = phase.tstart;
        b = a+phase.tdelta;
        if isa(a,'tomSym') || isa(b,'tomSym')
            if ~isempty(lastTomSymSolution)
                a = subs(a,lastTomSymSolution);
                b = subs(b,lastTomSymSolution);
            else
                error('EZPLOT could not find an ezsolve solution');
            end
        end
    else
        error('EZPLOT boundaries A and B are missing.');
    end
    nab = 0;
    coll = true;
    npts = max(100,3*length(phase.cpoints));
end

xx = linspace(a,b,npts);

if coll
    if ~isequal(s, phase.t)
        error('Free variable does not match the one used in phase.')
    end
    % Collocation points + intial and endpoint
    xxc = phase.tstart + ...
        phase.tdelta*unique([0; phase.cpoints(:); 1]);
    if isa(xxc,'tomSym')
        if ~isempty(lastTomSymSolution)
            xxc = subs(xxc,lastTomSymSolution);
        else
            error('EZPLOT could not find an ezsolve solution');
        end
    end
    xx  = unique([xx'; xxc]);
    xc  = atPoints(phase,xxc,x);
    if ny
        yc  = atPoints(phase,xxc,y);
    end
end

xp = zeros(length(xx),numel(x));

for k=1:length(xx)
    xp(k,:) = vec(subs(x,s,xx(k)))';
end
    
if ny
    yp = zeros(length(xx),numel(y));
    for k=1:length(xx)
        yp(k,:) = vec(subs(y,s,xx(k)))';
    end
    if coll
        h1 = plot(xc,yc,'*',varargin{ny+nab+2:end},xp,yp,'k-');
    else
        h1 = plot(xp,yp,varargin{ny+nab+2:end});
    end
else
    if coll
        h1 = plot(xxc,xc,'*',varargin{ny+nab+2:end},xx,xp,'k-');
    else
        h1 = plot(xx,xp,varargin{ny+nab+2:end});
    end
end

if nargout>=1
    % Only return the handle if it was asked for
    h = h1;
end
