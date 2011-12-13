% uc_hs_test.m
%
% function A = uc_hs_test(probnrs,algnrs)
%
% uc_hs_test runs the unconstrained Schittkowski testproblems defined in
% uc_hs_prob with different solvers and stores the number of iterations
% function and gradient evaluations.
%
% INPUT:
% probNumbs   Problems to be run (in vector), default: all
% algnrs      Algorithms be used (in vector), default: algnrs = [1 2 3 4 9 10 12]
%             alg = 0, ucSolve Newton
%             alg = 1, ucSolve BFGS
%             alg = 2, ucSolve Inv BFGS
%             alg = 3, ucSolve Inv DFP
%             alg = 4, ucSolve DFP
%             alg = 5, ucSolve Fletcher-Reeves CG
%             alg = 6, ucSolve Polak-Ribiere CG
%             alg = 7, Fletcher conjugate descent CG-method
%             alg = 8, sTrustR
%             alg = 9, CONSTR
%             alg = 10, MINOS
%             alg = 11, NPSOL
%
% OUTPUT:
% A           Matrix containing structs with info.

function A = uhs_test2(probNumbs,algnrs)

if nargin < 2
    algnrs = [];
    if nargin < 1
        probNumbs = [];
    end
end

probFile  = 'uhs_prob';
probnames = feval(probFile); % Names of available problems

if isempty(probNumbs), probNumbs = 1:size(probnames,1); end
if isempty(algnrs), algnrs = [1 2 3 4 9 10 12]; end

pNames = probnames;
pNames=str2mat(pNames,'Average');
pNames=str2mat(pNames,'Failures');
probnames = pNames;

PriLevOpt = 0;
wait      = 0;
PriLev    = 0;

% Result flags
FAIL    = 0;
MAXITER = 1;
OK      = 2;

for j = 1:length(algnrs)
    alg = algnrs(j);
    if alg == 0; % ucSolve Newton
        Result = TestRun('ucSolve',0,probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 1; % ucSolve BFGS
        Result = TestRun('ucSolve',1,probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 2; % ucSolve Inv BFGS
        Result = TestRun('ucSolve',2,probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 3; % ucSolve Inv DFP
        Result = TestRun('ucSolve',3,probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 4; % ucSolve DFP
        Result = TestRun('ucSolve',4,probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 5; % ucSolve F-R CG
        Result = TestRun('ucSolve',5,probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 6; % ucSolve P-R CG
        Result = TestRun('ucSolve',6,probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 7; % ucSolve F CG
        Result = TestRun('ucSolve',7,probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 8; % sTrustR
        Result = TestRun('sTrustR',[],probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 9; % CONSTR
        Result = TestRun('CONSTR',[],probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 10; % MINOS
        Result = TestRun('MINOS',[],probFile,probNumbs,PriLevOpt,wait,PriLev);
    elseif alg == 11; % NPSOL
        Result = TestRun('NPSOL',[],probFile,probNumbs,PriLevOpt,wait,PriLev);
    end

    its_sum = 0;
    n_f_sum = 0;
    n_g_sum = 0;
    count = 0;
    for i=1:length(probNumbs)
        fstar = Result(i).f_k;
        f_opt = Result(i).Prob.f_opt;
        n_f   = Result(i).FuncEv;
        n_g   = Result(i).GradEv;
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

            if its >= Result(i).Prob.optParam.MaxIter
                flag = MAXITER;
            else
                flag = OK;
            end
        end

        tempstruct = struct('f',{fstar},'x',{x_k},'n_f',{n_f},'n_g',{n_g},...
            'its',{its},'flag',{flag});
        A(i,j) = tempstruct;
    end
    if count > 0
        its_mean = round(its_sum/count);
        n_f_mean = round(n_f_sum/count);
        n_g_mean = round(n_g_sum/count);
    else
        its_mean = 0;
        n_f_mean = 0;
        n_g_mean = 0;
    end

    % Average values
    tempstruct = struct('f',{fstar},'x',{x_k},'n_f',{n_f_mean},'n_g',{n_g_mean},...
        'its',{its_mean},'flag',{OK});
    A(i+1,j) = tempstruct;

    % Failures
    failures = length(probNumbs)-count;
    tempstruct = struct('f',{fstar},'x',{x_k},'n_f',{failures},'n_g',{failures},...
        'its',{failures},'flag',{OK});
    A(i+2,j) = tempstruct;
end

% CREATE TABLES
algnames  = [{'ucSolve Newton'},{'ucSolve BFGS'},{'ucSolve Inv BFGS'},...
    {'ucSolve Inv DFP'},{'ucSolve DFP'},{'ucSolve F-R CG'},...
    {'ucSolve P-R CG'},{'ucSolve F CG'},{'sTrustR'},{'CONSTR'},...
    {'MINOS'},{'NPSOL'}];

for j = 1:length(algnrs)
    for i = 1:length(probNumbs)+2
        %B(i,j) = A(i,j).n_r;
        if A(i,j).flag == FAIL
            tempstring = '---';
            %HKH
            % tempstring = '\ast';
            % tempstring = '\star';
            % tempstring = '\dagger';
            % tempstring = '\ddagger';
            % tempstring = '\circ';
            % tempstring = '\bullet';
            % tempstring = '\diamond';
            % tempstring = '\odot';
            % tempstring = '\oslash';
            % tempstring = '\otimes';
            % tempstring = '\ominus';
            % tempstring = '\oplus';
        else
            if A(i,j).flag == MAXITER
                tsrt = '$\ast$';
            else
                tsrt = ' ';
            end
            if i==length(probNumbs)+2
                tempstring = sprintf('%4d %s',A(i,j).n_f,tsrt);
            else
                tempstring = sprintf('%3d%4d%4d% %s',A(i,j).its,A(i,j).n_f,A(i,j).n_g,tsrt);
            end
        end
        tempstruct = struct('text',{tempstring});
        B(i,j) = tempstruct;
    end
end
frame = 0; % No frame
Tsize = 4;
lrc   = 3;
finame  = 'uc_hs_comp';
decs = zeros(1,length(algnrs));
%caption = 'Comparison of algorithmic performance for constrained HS problems.';
caption = [];
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

% function maketabl(A, frame, Tsize, lrc, finame, decs, caption, head, ...
%                   subhead, lStr, rStr, Note, hcols)
%
% maketabl. Make tex file with floating table environment from matrix A
%           Label is set the same as the name of the file
%           The table is centered on the page (and floating)
%
% A       Matrix, m by n. If value NaN, nothing is displayed
%         If A is structure, A.text is displayed.
% frame   Type of table.
%         0 = No frame
%         1 = Frame around the table
%         2 = Also frame around the header
%         3 = All entries are framed
% Tsize   Type size of table. Tsize < 0 ==> Set no size, use LaTeX default.
%         0 = tiny              1 = scriptsize (default)    2 = footnotesize
%         3 = small             4 = normalsize              5 = large
% lrc     Adjustment for each column element. Vector is expanded to correct
%         length
%         = 1 left  adjusted
%         = 2 right adjusted
%         = 3 centered
% finame  Name of tex file (without extension .tex)
% decs    Output col j in A with decs(j) decimals. decs(j) < 0 gives integer.
%         Vector is expanded to correct length. Default is integer output.
% caption Table header text
% head    Header. If empty, no header is displayed
% subhead Cell array with sub-header strings. size(subhead,2)==size(A,2)
% lStr    If not empty, use this string matrix as first column, to the left
%         size(lStr,1)==size(subhead,1) + size(A,1))
% rStr    If not empty, use this string matrix as last column, to the right
%         size(rStr,1)==size(subhead,1) + size(A,1))
% Note    Three dimensional matrix (m,n,# of notes). If Note(i,j,k)=1,then
%         the number k is set as a raised note (in \tiny) for element(i,j)
% hcols   Number of columns to use. If hcols has length 2, Use
%         columns hcols(1)-hcols(2) for header. NOT IMPLEMENTED YET
%
% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: hkh@tomopt.com
% % Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.0.2$
% Written May 1, 1997. Last modified Apr 1, 1997.
%
% MB modified for cls_test where A is a matrix of structures containing
% a string. Mar 21 1998