% Routine to be called after optimization.
%
% function Result=endSolve(Prob,Result)
%
% Catch up search directions and line search steps (optionally).
% Computing CPU time and real time elapsed

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1997-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written June 1, 1997.   Last modified Sept 18, 2009.

function Result=endSolve(Prob,Result)

global SEP_z SEP_Jz

if isempty(Result)
    return
end
TIME0 = Prob.TIME0;

if ~isempty(TIME0)
    Result.CPUtime = cputime-TIME0;
else
    Result.CPUtime = NaN;
end
TIME1 = Prob.TIME1;
if ~isempty(TIME1)
    Result.REALtime = etime(clock,TIME1);
else
    Result.REALtime = NaN;
end

% --------
% endbuild
% --------
% Create search directions, step lengths and variable limits for
% solver routines. Call aVbuild and reset flag BUILDP.

global alphaV p_dx F_X BUILDP
if BUILDP == 1
    alphaV=aVbuild(alphaV,size(p_dx,2));
end
BUILDP=0;

% End of endbuild script

if isfield(Prob.LS,'SepAlg')
    if Prob.LS.SepAlg > 0
        Result.SepLS.z = SEP_z;
        Result.SepLS.Jz= SEP_Jz;
    end
end

% Save iterations in result struct
% global p_dx alphaV F_X
global X_min X_max

Result.p_dx=p_dx;
Result.alphaV=alphaV;
Result.x_min=X_min;
Result.x_max=X_max;
Result.F_X=F_X;

% Put Prob into Result struct if this is not already done
if ~isfield(Result,'Prob')
    Result.Prob=Prob;
end

global n_f n_g n_H
global n_c n_dc n_d2c
global n_r n_J

if isempty(Result.FuncEv) | Result.FuncEv==0
    Result.FuncEv=n_f;
end
if isempty(Result.ConstrEv) | Result.ConstrEv==0
    Result.ConstrEv=n_c;
end
if isempty(Result.ResEv) | Result.ResEv==0
    Result.ResEv=n_r;
end
if isempty(Result.GradEv) | Result.GradEv==0
    Result.GradEv=n_g;
end
if isempty(Result.HessEv) | Result.HessEv==0
    Result.HessEv=n_H;
end
if isempty(Result.ConJacEv) | Result.ConJacEv==0
    Result.ConJacEv=n_dc;
end
if isempty(Result.ConHessEv) | Result.ConHessEv==0
    Result.ConHessEv=n_d2c;
end
if isempty(Result.JacEv) | Result.JacEv==0
    Result.JacEv=n_J;
end

% Assign f_k, g_k, H_k and so on values.

% Avoid unnecessary computations for solvers of type:
% if Result.solvType ~= checkType('glb')  9
% if Result.solvType ~= checkType('glc') 10
% if Result.solvType ~= checkType('cgo') 15

if ~isempty(Result.x_k)
    x_k = Result.x_k(:,1);
    if Result.FuncEv > 0 & isempty(Result.f_k)
        Result.f_k = nlp_f(x_k, Prob);
    end
    if Result.ResEv > 0 & isempty(Result.r_k)
        Result.r_k = nlp_r(x_k, Prob);
    end
    if isempty(Result.c_k) & Result.ConstrEv > 0
       nxk = size(x_k,2);
       Result.c_k       = zeros(Prob.mNonLin,nxk);
       for i=1:nxk
           Result.c_k(:,i) = nlp_c(x_k(:,i),Prob);
       end
       % Result.c_k = nlp_c(x_k, Prob);
    end
    if ~any(Result.solvType == [9 10 11 12 15])
        if Result.GradEv > 0 & isempty(Result.g_k)
            Result.g_k = nlp_g(x_k, Prob);
        end
        if Result.HessEv > 0 & isempty(Result.H_k)
            Result.H_k = nlp_H(x_k, Prob);
        end
        if Result.ConJacEv > 0 & isempty(Result.cJac)
            Result.cJac = nlp_dc(x_k, Prob);
        end

        % Some special treatment with d2L:
        %  If constraints were used: make sure there are lagrange
        %  multipliers.
        % HKH - do not compute 2nd order information if empty H-function
        if isempty(Result.d2L) & ~isempty(Prob.FUNCS.H)
            if Result.ConHessEv > 0 & ~isempty(Result.v_k)
                idx = length(Result.v_k)-Prob.mNonLin+1;
                if idx > length(Result.v_k)
                    error(['No lagrange multipliers for non-linear constraints ' ...
                        'found']);
                end
                lam = Result.v_k(idx:end);
                Result.d2L = nlp_d2L(x_k, lam, 2, Prob);
            elseif Result.HessEv > 0
                Result.d2L = nlp_d2L(x_k, [], 2, Prob);
            end
        end
        if Result.JacEv > 0 & isempty(Result.J_k)
            Result.J_k = nlp_J(x_k, Prob);
        end
    end
end

% Add constant f term Prob.fConstant to f_k
Result.f_k = Result.f_k + Prob.fConstant;

% MODIFICATION LOG
%
% 020702 hkh  Return directly if Result is empty
% 040402 hkh  Use fields TIME0 and TIME1 in Prob, instead of globals
% 040407 hkh  Add code to set ConJacEv and ConHessEv
% 050223 frhe Added code to set f_k, g_k and so on
% 050414 hkh  Check if Result.x_k is []
% 050422 hkh  Add constant term Prob.fConstant to Result.f_k
% 080227 hkh  Avoid gradient/Hessian estimates for glb,glc,cgo types
% 080310 hkh  Do not compute 2nd order information if empty H-function
% 080607 hkh  Minor cleanup
% 090813 med  mlint check
% 090919 hkh  Avoid gradient/Hessian call for MIQP and MINLP (type 11 and 12)
% 090919 hkh  Use x_k(:,1), GO solvers might return several x vectors
