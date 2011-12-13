% con_hs_test.m
%
% function A = con_hs_test(probnrs,algnrs)
%
% con_hs_test runs the constrained Schittkowski testproblems defined in
% con_hs_prob with different solvers and stores the number of iterations
% function, gradient and constraint evaluations.
%
% INPUT:
% probNumbs   Problems to be run (in vector), default: all
% algnrs      Algorithms be used (in vector), default: algnrs = [0 1 3 4 7 8]
%             alg = 0, nlpSolve
%             alg = 1, nlpSolve BFGS
%             alg = 2, sTrustR
%             alg = 3, conSolve
%             alg = 4, conSolve BFGS
%             alg = 5, SNOPT
%             alg = 6, MINOS
%             alg = 7, NPSOL
%             alg = 8, FMINCON
%
% OUTPUT:
% A           Matrix containing structs with info.

function A = chs_test(probNumbs,algnrs)

if nargin < 2
    algnrs = [];
    if nargin < 1
        probNumbs = [];
    end
end

probFile  = 'chs_prob';
probnames = feval(probFile); % Names of available problems

if isempty(probNumbs), probNumbs = 1:size(probnames,1); end
if isempty(algnrs), algnrs = [0 1 3 4 7 8]; end

pNames = probnames;
pNames=str2mat(pNames,'Average');
pNames=str2mat(pNames,'Failures');
probnames = pNames;

PriLevOpt = 0;
wait      = 0;
PriLev    = 2;

% Result flags
FAIL    = 0;
MAXITER = 1;
OK      = 2;

for j = 1:length(algnrs)
    alg = algnrs(j);
    if alg == 0; % nlpSolve
        Result = runtest('nlpSolve',0,probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 1; % nlpSolve BFGS
        Result = runtest('nlpSolve',1,probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 2; % sTrustR
        Result = runtest('sTrustR',[],probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 3; % conSolve
        Result = runtest('conSolve',0,probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 4; % conSolve BFGS
        Result = runtest('conSolve',1,probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 5; % SNOPT
        Result = runtest('SNOPT',[],probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 6; % MINOS
        Result = runtest('MINOS',[],probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 7; % NPSL
        Result = runtest('NPSOL',[],probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 8; % FMINCON
        Result = runtest('FMINCON',[],probFile,probNumbs,PriLevOpt,wait,PriLev);
    end

    its_sum = 0;
    n_f_sum = 0;
    n_g_sum = 0;
    n_c_sum = 0;
    count = 0;
    for i=1:length(probNumbs)
        fstar = Result(i).f_k;
        f_opt = Result(i).Prob.f_opt;
        n_f   = Result(i).FuncEv;
        n_g   = Result(i).GradEv;
        n_c   = Result(i).ConstrEv;
        if isempty(Result(i).Iter)
            its = size(Result(i).p_dx,2); % OPTIM leastsq algs.
        else
            its   = Result(i).Iter;
        end
        x_k   = Result(i).x_k;


        if abs(fstar-f_opt)/max(1e-6,abs(f_opt)) > 1e-4 % Optimum function value not reached
            flag = FAIL;
        elseif 0 % Constraint satisfaction should be tested here
            flag = FAIL;
        else
            % If not failure compute average number of its, resev. Jacev.
            count = count + 1;
            its_sum = its_sum + its;
            n_f_sum = n_f_sum + n_f;
            n_g_sum = n_g_sum + n_g;
            n_c_sum = n_c_sum + n_c;

            if its >= Result(i).Prob.optParam.MaxIter
                flag = MAXITER;
            else
                flag = OK;
            end
        end

        tempstruct = struct('f',{fstar},'x',{x_k},'n_f',{n_f},'n_g',{n_g},...
            'n_c',{n_c},'its',{its},'flag',{flag});
        A(i,j) = tempstruct;
    end
    its_mean = round(its_sum/count);
    n_f_mean = round(n_f_sum/count);
    n_g_mean = round(n_g_sum/count);
    n_c_mean = round(n_c_sum/count);

    % Average values
    tempstruct = struct('f',{fstar},'x',{x_k},'n_f',{n_f_mean},'n_g',{n_g_mean},...
        'n_c',{n_c_mean},'its',{its_mean},'flag',{OK});
    A(i+1,j) = tempstruct;

    % Failures
    failures = length(probNumbs)-count;
    tempstruct = struct('f',{fstar},'x',{x_k},'n_f',{failures},'n_g',{failures},...
        'n_c',{failures},'its',{failures},'flag',{OK});
    A(i+2,j) = tempstruct;
end

% CREATE TABLES
algnames  = [{'nlpSolve'},{'nlpSolve BFGS'},{'sTrustR'},{'conSolve'},{'conSolve BFGS'},...
    {'SNOPT'},{'MINOS'},{'NPSOL'},{'FMINCON'}];

for j = 1:length(algnrs)
    for i = 1:length(probNumbs)+2
        if A(i,j).flag == FAIL
            tempstring = '---';
        else
            if A(i,j).flag == MAXITER
                tsrt = '$\ast$';
            else
                tsrt = ' ';
            end
            if i==length(probNumbs)+2
                tempstring = sprintf('%4d %s',A(i,j).n_f,tsrt);
            else
                tempstring = sprintf('%3d%4d%4d%4d %s',...
                    A(i,j).its,A(i,j).n_f,A(i,j).n_g,A(i,j).n_c,tsrt);
            end
        end
        tempstruct = struct('text',{tempstring});
        B(i,j) = tempstruct;
    end
end
frame = 2; % No frame
Tsize = 4;
lrc   = 3;
finame  = 'con_hs_comp';
decs = zeros(1,length(algnrs));
caption = 'Comparison of algorithmic performance for constrained HS problems.';
head = [];
subhead = algnames(algnrs+1);
dummystring =  ' ';
for i=1:min(10,size(probnames,2))-1
    dummystring = [dummystring,' '];
end
lStr = [dummystring;probnames([probNumbs size(probnames,1)-1 size(probnames,1)],1:min(10,size(probnames,2)))];
rStr = [];
Note = [];
hcols = [];

maketabl(B, frame, Tsize, lrc, finame, decs, caption, head, ...
    subhead, lStr, rStr, Note, hcols)