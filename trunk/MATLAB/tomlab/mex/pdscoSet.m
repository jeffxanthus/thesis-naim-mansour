function options = pdscoSet(varargin)
%pdscoSet creates or alters an options structure for pdsco.m.
%
%   options = pdscoSet('PARAM1',VALUE1,'PARAM2',VALUE2,...);
%   creates a structure with the specified named parameters and values.
%   Any unspecified parameters are set to [].
%   It is sufficient to type only the leading characters that uniquely
%   identify the parameter.  Case is ignored for parameter names.
%
%   NOTE: For values that are strings, correct case and the complete string
%   are required; if an invalid string is provided, the default is used.
%
%   options = pdscoSet(oldopts,'PARAM1',VALUE1,...);
%   creates a copy of oldopts with the named parameters reset to the
%   specified values.
%
%   options = pdscoSet(oldopts,newopts);
%   combines an existing structure oldopts with a new structure newopts.
%   Any parameters in newopts with non-empty values overwrite the
%   corresponding old parameters in oldopts.
%
%   pdscoSet with no input arguments and no output arguments displays all
%   parameter names and their possible values.
%
%   options = pdscoSet (with no input arguments) creates a structure
%   where all the fields are set to their default values.
%
% options.gamma       may be zero or positive.  A typical value is
%                     gamma = 1.0e-4.
% options.delta       must be positive.  Again, a typical value is
%                     delta = 1.0e-4.
%                     Values smaller than 1.0e-6 or 1.0e-7 may affect
%                     numerical reliability.  Increasing delta usually
%                     improves convergence.
%                     For least-squares applications, delta = 1.0 is
%                     appropriate.
% options.MaxIter     is the maximum iterations of the primal-dual
%                     barrier method.
% options.FeaTol      is the accuracy for satisfying Ax + r = b and
%                     A'y + z = g.
% options.OptTol      is the accuracy for satisfying x.*z = 0.
% options.StepTol     is how close each step in x and z may be to
%                     reaching a bound.
% options.x0min       Minimum size of x0 relative to max(x0).
% options.z0min       Minimum size of z0 relative to max(z0).
% options.mu0         Scales initial mu.
% options.Method      1=Cholesky    2=QR    3=LSQR or SYMMLQ
%   4 = Tomlab Tlsqr, special PDSCO interface avoiding any callbacks (default)
% options.LSproblem   LS problem for each search direction.
% options.LSQRMaxIter * min(m,n) is the maximum LSQR (cg) iterations.
% options.LSQRatol1   is the starting value of the LSQR accuracy
%                     tolerance "atol".
% options.LSQRatol2   is the smallest value atol is reduced to.
% options.LSQRconlim  shuts LSQR down early if its matrix is ill-conditioned.
% options.wait = 0    means solve the problem with default internal parameters;
%              = 1    means pause to allow interactive resetting of parameters.


% pdscoSet.m is derived from optimset.m (Revision 1.14, 1998/08/17)
% in the Optimization Toolbox of The MathWorks, Inc.
%
% 28 Sep 2000: First version of pdscoSet.m.
%              Michael Saunders, SOL, Stanford University.
% 25 Jan 2003: More LSmethod options added for use in Tomlab (K. Holmstrom)
% 30 Aug 2004: LSmethod changed to Method (A. Goran/K. Holmstrom)
% 01 Aug 2005: isstr changed to ischar

if (nargin == 0)        % Set default options.
    defoptions.gamma        =  1e-4;
    defoptions.delta        =  1e-4;
    defoptions.MaxIter      =    30;
    defoptions.FeaTol       =  1e-6;
    defoptions.OptTol       =  1e-6;
    defoptions.StepTol      =  0.99;
    defoptions.x0min        =   1.0;  % 1.0 for cold starts?
    defoptions.z0min        =   1.0;  % 0.1 for warm starts?
    defoptions.mu0          =  1e-4;  % << 1.0!
    defoptions.Method       =     4;
    defoptions.LSproblem    =     1;
    defoptions.LSQRMaxIter  =  10.0;
    defoptions.LSQRatol1    = 1e-04;
    defoptions.LSQRatol2    = 1e-15;  % Not used
    defoptions.LSQRconlim   = 1e+12;  % Somewhere between e+8 and e+16
    defoptions.wait         =     0;
    defoptions.NOTE         = 'LSQRMaxIter is scaled by the matrix dimension';

    if (nargout == 0)    % Display options.
       disp('pdsco default options:')
       disp( defoptions )
    else
       options = defoptions;
    end
    return;
end

Names = ...
[
    'gamma      '
    'delta      '
    'MaxIter    '
    'FeaTol     '
    'OptTol     '
    'StepTol    '
    'x0min      '
    'z0min      '
    'mu0        '
    'Method     '
    'LSproblem  '
    'LSQRMaxIter'
    'LSQRatol1  '
    'LSQRatol2  '
    'LSQRconlim '
    'wait       '
    'NOTE       '
];
m     = size (Names,1);
names = lower(Names);

% The remaining clever stuff is from optimset.m.

% Combine all leading options structures o1, o2, ... in pdscoSet(o1,o2,...).
options = [];
for j = 1:m
    eval(['options.' Names(j,:) '= [];']);
end
i = 1;
while i <= nargin
    arg = varargin{i};
    if ischar(arg)                         % arg is an option name
       break;
    end
    if ~isempty(arg)                      % [] is a valid options argument
       if ~isa(arg,'struct')
          error(sprintf(['Expected argument %d to be a ' ...
                'string parameter name ' ...
                'or an options structure\ncreated with pdscoSet.'], i));
       end
       for j = 1:m
          if any(strcmp(fieldnames(arg),deblank(Names(j,:))))
             eval(['val = arg.' Names(j,:) ';']);
          else
             val = [];
          end
          if ~isempty(val)
             eval(['options.' Names(j,:) '= val;']);
          end
       end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end

expectval = 0;                          % start expecting a name, not a value

while i <= nargin
    arg = varargin{i};

    if ~expectval
       if ~ischar(arg)
          error(sprintf(['Expected argument %d to be a ' ...
                         'string parameter name.'], i));
       end

       lowArg = lower(arg);
       j = strmatch(lowArg,names);
       if isempty(j)                       % if no matches
          error(sprintf('Unrecognized parameter name ''%s''.', arg));
       elseif length(j) > 1                % if more than one match
          % Check for any exact matches
          % (in case any names are subsets of others)
          k = strmatch(lowArg,names,'exact');
          if length(k) == 1
             j = k;
          else
             msg = sprintf('Ambiguous parameter name ''%s'' ', arg);
             msg = [msg '(' deblank(Names(j(1),:))];
             for k = j(2:length(j))'
                msg = [msg ', ' deblank(Names(k,:))];
             end
             msg = sprintf('%s).', msg);
             error(msg);
          end
       end
    else
       eval(['options.' Names(j,:) '= arg;']);
    end

    expectval = ~expectval;
    i = i + 1;
end

if expectval
    error(sprintf('Expected value for parameter ''%s''.', arg));
end
