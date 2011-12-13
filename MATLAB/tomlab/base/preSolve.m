% preSolve.m
%
% function Prob = preSolve(Prob)
%
% preSolve is a presolve procedure for linear constraints and simple bounds.
%
% preSolve calls the subroutines:
%
% clean        Calls the routines emptyrow.m, r_rw_sng.m and el_cnsts.m	who
%              are executed until no changes in constraints or bounds are made.
%
% emptyrow     Eliminates empty rows from the constraint matrix.
%
% r_rw_sng     Removes row singletons.
%
% el_cnsts     Improves variable bounds and uses them to eliminate
%              redundant and forcing constraints.
%
% mksp         Search for sparsity patterns and eliminates some sort
%              of linear dependent rows.
%
% INPUT PARAMETERS
% Prob    Structure, where the following variables are used:
%         A         Linear constraint matrix
%         b_L       Lower bounds for linear constraints
%         b_U       Upper bounds for linear constraints
%         x_L       Lower bounds for x
%         x_U       Upper bounds for x
%         PriLevOpt Printing level
% ---------------------------------------
% MIP         Structure in Prob, Prob.MIP
% ---------------------------------------
%             Defines integer optimization parameters. Fields used:
%   IntVars:  
%             If empty, all variables are assumed non-integer 
%             If islogical(IntVars) (=all elements are 0/1), then
%             1 = integer variable, 0 = continuous variable.
%             If any element >1, IntVars is the indices for integer variables
%
% OUTPUT PARAMETERS
% Prob    Structure, where the following variables are changed:
%         A       Linear constraint matrix
%         b_L     Lower bounds for linear constraints, -inf if made redundant
%         b_U     Upper bounds for linear constraints,  inf if made redundant
%         If both b_L and b_U are infinite, the constraint is eliminated.
%         The original matrix A and b_L and b_U are stored in Prob.preSolve
%         x_L     New lower bounds for x
%         x_U     New upper bounds for x
%         
%
% New output fields in structure Prob.preSolve:
%         idxRed  Indices for the redundant linear constraints
%                 These constraints now have b_L = -inf, b_U = inf
%         idxL    The lower bounds made redundant by preSolve
%         idxU    The upper bounds made redundant by preSolve
%         A       Linear constraint matrix (Original)
%         b_L     Lower bounds for linear constraints (Original)
%         b_U     Upper bounds for linear constraints (Original)
%         x_L     Lower bounds for x (Original)
%         x_U     Upper bounds for x (Original)

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2007 by Tomlab Optimization Inc., $Release: 5.8.0$
% Written Januari 25, 1998.  Last modified Feb 22, 2007.


%****************************** PRESOLVE ******************************

function Prob = preSolve(Prob)

global feasible
global aeps
global beps

[m,n]  = size(Prob.A);

if m==0
   return
end
% Integer variables
IntVars  = DefPar(Prob.MIP,'IntVars',[]);

% Logical vector for integers
IV = false(n,1);

if isempty(IntVars)
   % No binary variables B or integer variables of type I
   MIP = 0;
elseif any(IntVars==0) | all(IntVars==1)
   % Assume binary logical vector given
   IV(1:length(IntVars)) = logical(IntVars);
   MIP = 1;
else
   if any(IntVars < 1 | IntVars > n)
      error('presolve: Illegal IntVars vector');
   end
   IV(IntVars)=1;
   MIP = 1;
end

% Pick up input variables from Prob structure
%B   = Prob.A;
b_L = Prob.b_L(:);
b_U = Prob.b_U(:);
x_L = Prob.x_L(:);
x_U = Prob.x_U(:);
if MIP 
   x_L(IV) = round(x_L(IV)+0.49999);
   x_U(IV) = round(x_U(IV)-0.49999);
end
if issparse(Prob.A)
   SPARSE=1;
   nnzA = nnz(Prob.A);
else
   SPARSE=0;
end

PriLev = Prob.PriLevOpt;       % Print level

if isempty(x_U), 
   x_U = inf*ones(n,1); 
   Prob.x_U = x_U;
end
if isempty(x_L) 
   x_L =-inf*ones(n,1); 
   Prob.x_L = x_L;
end
if isempty(b_U) 
   b_U = inf*ones(m,1); 
   Prob.b_U = b_U;
end
if isempty(b_L) 
   b_L =-inf*ones(m,1); 
   Prob.b_L = b_L;
end
Prob.preSolve.x_L = x_L;
Prob.preSolve.x_U = x_U;
Prob.preSolve.A   = Prob.A;
Prob.preSolve.b_L = b_L;
Prob.preSolve.b_U = b_U;

feasible = 1;                     % Feasibility flag
beps = eps*max(max(abs(Prob.A))); % Right hand side tolerance
aeps = eps;                       % Matrix entry tolerance

% Put equality constraints first in A and turn >= inequalities
% to <= inequalities.
Eq = find(b_L==b_U & ~isinf(b_L));
Lq = find(~isinf(b_L) & (b_L~=b_U));
Uq = find(~isinf(b_U) & (b_L~=b_U));
%Lq = find((b_L>-inf) & (b_L~=b_U));
%Uq = find((b_U< inf) & (b_L~=b_U));
if SPARSE
   A  = sparse([Prob.A(Eq,:);-Prob.A(Lq,:);Prob.A(Uq,:)]);
else
   A  = [Prob.A(Eq,:);-Prob.A(Lq,:);Prob.A(Uq,:)];
end

b  = [b_L(Eq);-b_L(Lq);b_U(Uq)];
EqCon = [ones(length(Eq),1);zeros(length(Lq)+length(Uq),1)];
   

% Call subroutines
[A, b, EqCon, x_L, x_U] = clean(A, b, EqCon, x_L, x_U, IV, PriLev);
[A, b, EqCon, mkspflag] = mksp(A, b, EqCon, IV);
if mkspflag
   [A, b, EqCon, x_L, x_U] = clean(A, b, EqCon, x_L, x_U, IV, PriLev);
end
   
% Reconstruct constraint matrix and bounds.   
if SPARSE
   % A_new   = spalloc(m,n,nnz(A)sum(A~=0));
   %A_new   = spalloc(m,n,nnz(A));
   Prob.A  = spalloc(m,n,nnzA);
else
   Prob.A  = zeros(m,n);
end

% Rows corresponding to equality constraints
Prob.A(Eq,:) = A(1:length(Eq),:);
b_L(Eq)      = b(1:length(Eq));
b_U(Eq)      = b(1:length(Eq));

% Rows corresponding to >= inequality constraints
Prob.A(Lq,:)  = -A(length(Eq)+1:length(Eq)+length(Lq),:);
b_L(Lq)       = -b(length(Eq)+1:length(Eq)+length(Lq));

% Rows corresponding to <= inequality constraints
Prob.A(Uq,:)  = A(length(Eq)+length(Lq)+1:length(b),:);
b_U(Uq)       = b(length(Eq)+length(Lq)+1:length(b));

%idx1 = find(isnan(b_L));
%idx2 = find(isnan(b_U));
%b_L(idx2) = NaN;
%b_U(idx1) = NaN;

idxL = find(isnan(b_L));
idxU = find(isnan(b_U));
b_L(idxL)  = -inf;
b_U(idxU)  = inf;

% Redundant constraints have both bounds infinity
idxRed               = find(isinf(b_L) & isinf(b_U));
Prob.preSolve.idxRed = idxRed;
Prob.preSolve.idxL   = idxL;
Prob.preSolve.idxU   = idxU;

if PriLev > 0
   fprintf('\n');
   fprintf('preSolve Analysis ');
   fprintf('\n');
   fprintf('Variables   %d. ',n);
   fprintf('Constraints %d. --- ',m);
   fprintf('Redundant constraints %d ',length(idxRed));
   fprintf('\n');
   if PriLev > 1
      Prob.b_L
      b_L
   end

   ibL = find(Prob.b_L(:)~=b_L & ~isinf(b_L));
   ibU = find(Prob.b_U(:)~=b_U & ~isinf(b_U));
   if ~isempty(idxL)
      fprintf('Eliminated lower bounds for %d constraints\n',length(idxL));
   end
   if ~isempty(ibL)
      fprintf('Changed    lower bounds for %d constraints\n',length(ibL));
   end
   if ~isempty(idxU)
      fprintf('Eliminated upper bounds for %d constraints\n',length(idxU));
   end
   if ~isempty(ibU)
      fprintf('Changed    upper bounds for %d constraints\n',length(ibU));
   end
   fprintf('\n');
   ixL = find(Prob.x_L(:)~=x_L);
   ixU = find(Prob.x_U(:)~=x_U);
   if ~isempty(ixL)
      fprintf('Changed lower bounds for %d variables\n',length(ixL));
   end
   if ~isempty(ixU)
      fprintf('Changed upper bounds for %d variables\n',length(ixU));
   end
   fprintf('\n');
   if PriLev > 1
      if ~isempty(idxRed)
         xprinti(idxRed,'Redundant:');
      end
      if ~isempty(idxL)
         xprinti(idxL,'Low Elimin:');
      end
      if ~isempty(ibL)
         xprinti(ibL,'Low Change:');
         xprint(Prob.b_L(ibL),'Old:');
         xprint(b_L(ibL),'New:');
      end
      if ~isempty(idxU)
         xprinti(idxU,'Upp Elimin:');
      end
      if ~isempty(ibU)
         xprinti(ibU,'Upp Change:');
         xprint(Prob.b_U(ibU),'Old:');
         xprint(b_U(ibU),'New:');
      end
      if ~isempty(ixL)
         xprinti(ixL,'Lower:');
         xprint(Prob.x_L(ixL),'Old:');
         xprint(x_L(ixL),'New:');
      end
      if ~isempty(ixU)
         xprinti(ixU,'Upper:');
         xprint(Prob.x_U(ixU),'Old:');
         xprint(x_U(ixU),'New:');
      end
   end
end

%Prob.A = A_new;
%Prob.b_L = b_L;
%Prob.b_U = b_U;
Prob.x_L = x_L;
Prob.x_U = x_U;
if ~isempty(idxRed)
   ix         = ones(m,1);
   ix(idxRed) = 0;
   ixA        = find(ix);
   Prob.A     = Prob.A(ixA,:);
   Prob.b_L   = b_L(ixA);
   Prob.b_U   = b_U(ixA);
   Prob.mLin  = size(Prob.A,1);
   if PriLev > 0
      fprintf('Eliminated %d linear constraints\n',length(idxRed));
   end
end

%******************************* CLEAN ******************************
function [A, b, EqCon, x_L, x_U] = clean(A, b, EqCon, x_L, x_U, IV, PriLev)

global feasible

change = 1; % Flag if changes are made
while change & feasible 
   % Call emptyrow, r_rw_sng and el_cnsts until no changes are made
   change = 0;
   [A, b, EqCon] = emptyrow(A, b, EqCon, IV, PriLev);
   [A, b, EqCon, x_L, x_U, change] = r_rw_sng(A, b, EqCon,...
    x_L, x_U, change, IV, PriLev);
   [A, b, EqCon] = emptyrow(A, b, EqCon, IV, PriLev);
   [A, b, EqCon, x_L, x_U, change] = el_cnsts(A, b, EqCon,...
    x_L, x_U, change, IV, PriLev);
end



%****************************** EMPTYROW ******************************
function [A, b, EqCon] = emptyrow(A, b, EqCon, IV, PriLev)

global feasible
global aeps
global beps

if feasible
   [m,n]  = size(A);
   removerow = 0;      % Flag if any row is to be removed
   rowidx = ones(m,1); % Logic vector for constraints
  
   for i = 1:m
      if ~isnan(b(i)) & feasible
         if ~any(abs(A(i,:)) > aeps) % Empty row
            if (~EqCon(i) & b(i)<-beps) | ...
               (EqCon(i) & abs(b(i))>beps) % Infeasible
            %if ((EqCon(i)~=11) & b(i)<-beps) | ...
            %   ((EqCon(i)==11) & abs(b(i))>beps) % Infeasible
               feasible = 0;
               if PriLev>0 
                  fprintf('\n Infeasible due to constraint %d',i);
               end
            else % Redundant constraint
               rowidx(i) = 0;
               if PriLev > 4, fprintf('\n Row %d is empty',i); end
               removerow = 1;
            end
         end
      end
   end
   if removerow & feasible % Delete redundant constraints
      b(find(~rowidx))=NaN; 
   end
end



%****************************** R_RW_SNG ******************************
function [A, b, EqCon, x_L, x_U, change] = r_rw_sng(A, b, EqCon,...
          x_L, x_U, change, IV, PriLev)

global feasible    
global aeps
global beps

if feasible
   [m,n]  = size(A);
   removerow = 0; % Flag if any row is to be removed
   rowidx = ones(m,1); % Logic vector for constraints
   % HKH, add 3 lines
   AA = abs( A ) > aeps;
   AAsum = sum( AA' );
   %v  = ~isnan(b);

   for i = 1:m
   bi = b(i);
   if ~isnan(bi) & feasible
      %if (length( find( abs(A(i,:)) > aeps ) ) == 1) % Singleton row
      %   j = find( abs(A(i,:)) > aeps );
      if (AAsum(i) == 1) % Singleton row
         j = find( AA(i,:));
         %if EqCon(i)~=11 % Inequality constraint (<=)
         if ~EqCon(i) % Inequality constraint (<=)
            bound = bi/A(i,j);
            if A(i,j) > 0 % Bound is a upper bound
               if bound < x_L(j) % Infeasible
                  feasible = 0;
                  if PriLev > 0
                     fprintf('\n Infeasible due to constraint %d',i);
                     fprintf('\n Constraints gives upper bound ');
                     fprintf('below x_L for variable %d',j);
                  end
               else % Check if x_U is to be tightened
                  if x_U(j) > bound
                     x_U(j) = bound;
                     if IV(j), x_U(j) = round(x_U(j) - 0.49999); end
                     change = 1;
                  end                    
               end
            else % Bound is a lower bound
               if bound > x_U(j) % Infeasible
                  feasible = 0;
                  if PriLev > 0
                     fprintf('\n Infeasible due to constraint %d',i);
                     fprintf('\n Constraints gives lower bound ');
                     fprintf('above x_U for variable %d',j);
                  end
               else % Check if x_L is to be tightened
                  if x_L(j) < bound
                     x_L(j) = bound;
                     if IV(j), x_L(j) = round(x_L(j) + 0.49999); end
                     change = 1;
                  end
               end
            end
         else % Equality constraint
            fixvalue = bi/A(i,j);
            if ~( (x_L(j) <= fixvalue) & (fixvalue<=x_U(j)) ) % Infeasible
               feasible = 0;
               if PriLev > 0
                  fprintf('\n Infeasible due to constraint %d',i);
                  fprintf('\n Constraints fix variable %d beyond the bounds',j);
               end
            else % Set lower and upper bounds to fixvalue and eliminate
                 % variable j in other equations.
               x_L(j) = fixvalue; x_U(j) = fixvalue;
               % HKH Could put check on integer value here !!!
               % New code
               v = AA(:,j);
               v(i) = 0;
               ix = find(v);
               if ~isempty(ix)
                  b(ix) = b(ix) - A(ix,j)*fixvalue;
                  A(ix,j) = 0;
                  % Update logical matrix and sum of nonzeros
                  AA(ix,j) = 0;
                  AAsum(ix) = AAsum(ix)-1;
               end
               % End new code
               %for row = 1:m
               %   if (row ~= i) & (abs(A(row,j))>aeps)
               %      b(row) = b(row) - A(row,j)*fixvalue;
               %      A(row,j) = 0;
               %   end
               %end
            end
         end
         rowidx(i) = 0;
         removerow = 1;
         if PriLev > 4 
            fprintf('\n Row %d is singleton in variable %d',i,j);
         end
      end
   end
   end 
   if removerow & feasible % Delete redundant constraints
      b(find(~rowidx))=NaN; 
      change = 1;      
   end
end



%****************************** EL_CNSTS ******************************
function [A, b, EqCon, x_L, x_U, change] = el_cnsts(A, b, EqCon,...
          x_L, x_U, change, IV, PriLev)

global feasible
global aeps
global beps

if feasible   
   isxL = isfinite(x_L)';
   isxU = isfinite(x_U)';
   [m,n]  = size(A);
   removerow = 0; % Flag if any row is to be removed
   rowidx = ones(m,1); % Logic vector for constraints
   AA =  A > aeps;
   AN =  A < -aeps;
   for i = 1:m
   bi = b(i);
   if ~isnan(bi) & feasible
      if ~EqCon(i)  % Inequalities only this far
      %if (EqCon(i)~=11)  % Inequalities only this far
         blow = 0; % Lower bound for b(i)
         bup  = 0; % Upper bound for b(i)
         ix = find(AA(i,:));
         for k =1:length(ix)
             j = ix(k);
             aij = A(i,j);
             blow = blow + aij*x_L(j);
             bup  = bup  + aij*x_U(j);
         end
         ix = find(AN(i,:));
         for k =1:length(ix)
             j = ix(k);
             aij = A(i,j);
             blow = blow + aij*x_U(j);
             bup  = bup  + aij*x_L(j);
         end
         %for j = 1:n
         %    aij = A(i,j);
         %    if aij > aeps
         %       blow = blow + aij*x_L(j);
         %       bup  = bup  + aij*x_U(j);
         %    elseif aij < -aeps
         %       blow = blow + aij*x_U(j);
         %       bup  = bup  + aij*x_L(j);
         %    end
         %end

         % Four possibilities may occur
         if bi < blow % Infeasible
            feasible = 0;
            if PriLev > 0 
               fprintf('\n Infeasible due to constraint %d',i);
               fprintf('\n Variable bounds gives lower bound above b(%d)',i);
            end
         elseif abs( bi-blow ) < beps % blow equals b(i)
            for k =1:length(ix)
                j = ix(k);
                x_U(j) = x_L(j);
            end
            ix = find(AN(i,:));
            for k =1:length(ix)
                j = ix(k);
                x_L(j) = x_U(j);
            end
            %for j=1:n
            %    aij = A(i,j);
            %    if aij > aeps % Fix vars with A(i,j)>0 on lower bounds
            %       x_U(j) = x_L(j);
            %    elseif aij < -aeps % Fix vars with A(i,j)<0 on upper bounds
            %       x_L(j) = x_U(j);
            %    end
            %end
            rowidx(i) = 0;
            removerow = 1;
            if PriLev > 4 
               fprintf('\n Constraint %d is eliminated as a forcing one',i);
            end
         elseif bup <= bi % Redundant constraint
            rowidx(i) = 0;
            removerow = 1;
            if PriLev > 4 
               fprintf('\n Constraint %d is eliminated as a redundant one',i);
            end
         elseif ( blow-bi < -beps ) & ( bup-bi > beps ) 
            % Tightening variable bounds
            if isfinite( blow )
               ix = find(AA(i,:));
               for k =1:length(ix)
                   j = ix(k);
                   if isxL(j)
                      Uprime = x_L(j) + ( bi-blow )/A(i,j);
                      if Uprime < x_L(j) % Infeasible
                         feasible = 0;
                         if PriLev > 0 
                            fprintf('\n Infeasible due to constraint %d',i);
                            fprintf('\n New upper bound is below ');
                            fprintf('lower bound for variable %d',j);
                         end
                      elseif Uprime < x_U(j) % Constructive, tighten bound
                         x_U(j) = Uprime;
                         if IV(j), x_U(j) = round(x_U(j) - 0.49999); end
                         change = 1;
                      end  % Else new bound redundant
                   end
               end
               ix = find(AN(i,:));
               for k =1:length(ix)
                   j = ix(k);
                   if isxU(j)
                      Lprime = x_U(j) + ( bi-blow )/ A(i,j);
                      if Lprime > x_U(j) % Infeasible
                         feasible = 0;
                         if PriLev > 0 
                            fprintf('\n Infeasible due to constraint %d',i);
                            fprintf('\n New lower bound is above ');
                            fprintf('upper bound for variable %d',j);
                         end
                      elseif Lprime > x_L(j) % Constructive, tighten bound
                         x_L(j) = Lprime;
                         if IV(j), x_L(j) = round(x_L(j) + 0.49999); end
                         change = 1;
                      end % Else new bound redundant
                   end
               end
               %if 0
               %for j = 1:n
               %    aij = A(i,j);
               %    %if ( aij > aeps ) & isfinite(x_L(j)),zz=1; 
               %    %elseif ( aij < -aeps ) & isfinite(x_U(j))
               %    %        zz=1;
               %    %end
               %    if ( aij > aeps ) & isxL(j)
               %       Uprime = x_L(j) + ( bi-blow )/aij;
               %       if Uprime < x_L(j) % Infeasible
               %          feasible = 0;
               %          if PriLev > 0 
               %             fprintf('\n Infeasible due to constraint %d',i);
               %             fprintf('\n New upper bound is below ');
               %             fprintf('lower bound for variable %d',j);
               %          end
               %       elseif Uprime < x_U(j) % Constructive, tighten bound
               %          x_U(j) = Uprime;
               %          change = 1;
               %       end  % Else new bound redundant
               %    elseif ( aij < -aeps ) & isxU(j)
               %       Lprime = x_U(j) + ( bi-blow )/ aij;
               %       if Lprime > x_U(j) % Infeasible
               %          feasible = 0;
               %          if PriLev > 0 
               %             fprintf('\n Infeasible due to constraint %d',i);
               %             fprintf('\n New lower bound is above ');
               %             fprintf('upper bound for variable %d',j);
               %          end
               %       elseif Lprime > x_L(j) % Constructive, tighten bound
               %          x_L(j) = Lprime;
               %          change = 1;
               %       end % Else new bound redundant
               %    end
               %end
               %end
            end
            
            % noilb - Number Of Infinite Lower Bounds for variables 
            %         with positive entry in row i 
            % noiub - Number Of Infinite Upper Bounds for variables 
            %         with negative entry in row i 
            %noilb = sum(~isfinite(x_L(find( A(i,:) > aeps ))));
            %noiub = sum(~isfinite(x_U(find( A(i,:) <-aeps ))));
            %noilb2 = sum(~isxL(find( A(i,:) > aeps )));
            %noiub2 = sum(~isxU(find( A(i,:) <-aeps )));
            noilb = sum(~isxL & AA(i,:));
            noiub = sum(~isxU & AN(i,:));
            %if any(noilb2~=noilb), keyboard, end
            %if any(noiub2~=noiub), keyboard, end
            if  (noilb+noiub) == 1
               % Only one infinite bound for vars appearing in row i   
               if noilb == 1 % Infinite bound is a lower bound
                  %k = find(~isfinite(x_L(find( A(i,:) > aeps ))));
                  %k = find(~isxL(find( A(i,:) > aeps )));
                  k = find(~isxL & AA(i,:));
               else % Infinite bound is upper bound
                  %k = find(~isfinite(x_U(find( A(i,:) <-aeps ))));
                  %k = find(~isxU(find( A(i,:) <-aeps )));
                  k = find(~isxU & AN(i,:));
               end
               sum1=0;
               sum2=0;
               aik = A(i,k);
               if aik > 0 % New upper bound Uprime is to be determined
                  for j = 1:n
                      aij = A(i,j);
                      if aij > 0
                         if j ~= k
                            sum1 = sum1 + aij*x_L(j);
                         end
                      else
                         sum2 = sum2 + aij*x_U(j);
                      end
                  end
                  Uprime = ( bi -sum1 -sum2 )/aik;
                  if Uprime < x_L(j) % Infeasible
                     feasible = 0;
                     if PriLev > 0 
                        fprintf('\n Infeasible due to constraint %d',i);
                        fprintf('\n New upper bound is below ');
                        fprintf('lower bound for variable %d',j);
                     end
                  elseif Uprime < x_U(j) % Constructive, tighten bound
                     x_U(j) = Uprime;
                     if IV(j), x_U(j) = round(x_U(j) - 0.49999); end
                     change = 1;
                  end  % Else new bound redundant                  
               elseif aik < 0 % New lower bound Lprime is to be determined
               %else % New lower bound Lprime is to be determined
                  for j = 1:n
                     if aik > 0
                        sum1 = sum1 + A(i,j)*x_L(j);
                     else
                        if j ~= k
                           sum2 = sum2 + A(i,j)*x_U(j);
                        end
                     end
                  end
                  Lprime = ( bi -sum1 -sum2 )/aik;
                  if Lprime > x_U(j) % Infeasible
                     feasible = 0;
                     if PriLev > 0 
                        fprintf('\n Infeasible due to constraint %d',i);
                        fprintf('\n New lower bound is above ');
                        fprintf('upper bound for variable %d',j);
                     end
                  elseif Lprime > x_L(j) % Constructive, tighten bound
                     x_L(j) = Lprime;
                     if IV(j), x_L(j) = round(x_L(j) + 0.49999); end
                     change = 1;
                  end % Else new bound redundant                  
               end   
            end      
         end
      end
   end
   end
   if removerow & feasible % Delete redundant constraints
      b(find(~rowidx))=NaN; 
      change = 1;
   end   
end



%****************************** MKSP ******************************
function [A,b,EqCon,mkspflag] = mksp(A, b, EqCon, IV)

global feasible
global aeps
global beps

mkspflag = 0;	% Flag if A is made sparser
if feasible & ( length(EqCon) > 0 )  
   [m,n]  = size(A);
   
   AA = abs( A ) > aeps;
   v  = ~isnan(b);
   AAsum = sum( AA' );
   AAsum = AAsum(:);
   % candrowidx is a list of all row pivot candidates (equality type rows)
   %candrowidx = find( (EqCon==11) & (~isnan(b)) ); 
   %candRow = (EqCon==11) & v; 
   candRow = EqCon & v; 

   candrowidx = find( candRow ); 
   
   while length(candrowidx) > 0
      % candrowsum = number of nonzero entries in in the pivot candidate rows
      %candrowsum2 = sum( [abs( A(candrowidx,:) ) > aeps]' );
      %candrowsum = sum( AA(candrowidx,:)' );
      candrowsum = AAsum(candrowidx);
      temp = find(candrowsum > 0);
      %temp2 = find(candrowsum2 > 0);
      %if any(full(candrowsum2)~=full(candrowsum)), keyboard, end
      %if any(temp ~=temp2), keyboard, end
      if isempty(temp)
         break;
      end
      % Use row candrowidx(piv) as pivot row
      [dummy,pivtemp] = min(candrowsum(temp)); 
      piv = temp(pivtemp);      
      pivrow = candrowidx(piv);
      % Find all columns having an nonzero entry in pivrow
      %colidx = find(abs(A(pivrow,:)) > aeps);
      colidx = find(AA(pivrow,:));

      % Find pivcol, the shortest (least number of nonzero entries) 
      % column in colidx
      %tmpidx = find(~isnan(b)); 
      %colsum = sum( [abs( A(tmpidx,colidx) ) > aeps] );
      colsum = sum( AA(find(v),colidx) );
      [dummy,col] = min(colsum);
      pivcol = colidx(col);
      % Use A(pivrow,pivcol) as pivot element
      
      z = AA(pivrow,:);
      v(pivrow) = 0;
      %size(v)
      %size(AAsum)
      ix = find(v & AAsum >= AAsum(pivrow));
      v(pivrow) = 1;
      % Find those rows who are supersets of row pivrow
      %for i = 1:m
      %if ~isnan(b(i)) 
      %   if i ~= pivrow
      for j = 1:length(ix)
          i = ix(j);
            %if sum( (abs(A(pivrow,:))>aeps)<=(abs(A(i,:))>aeps) ) == n
            %z1 = sum( (abs(A(pivrow,:))>aeps)<=(abs(A(i,:))>aeps) ) == n;
            %z3 = all( (abs(A(pivrow,:))>aeps)<=(abs(A(i,:))>aeps) );
            %z2 = all( z <= AA(i,:) );
            %z4 = ~any( z > AA(i,:) );
            %if all( z <= AA(i,:) )
            if ~any( z > AA(i,:) )
               % Row i is a superset of piwrow
               % Eliminate entry pivcol in row i
               pivvalue = -A(i,pivcol)/A(pivrow,pivcol);
               A(i,:) = A(i,:) + pivvalue*A(pivrow,:);
               b(i) = b(i) + pivvalue*b(pivrow);
               % Reenter i to candrowidx if not already a member
               %if EqCon(i) == 11 
               if EqCon(i)
                  candRow(i) = 1;
                  %if ~any(candrowidx==i)
                  %   candrowidx = [candrowidx;i]; 
                  %end
               end
               mkspflag = 1;
               % HKH: Update logical matrix AA and sum of nonzeros for row i
               AA(i,:) = abs( A(i,:) ) > aeps;
               AAsum(i) = sum( AA(i,:) );
            end
         %end
      %end
      end 
      % Remove pivrow from candrowidx
      %candrowidx = [candrowidx(1:piv-1);candrowidx(piv+1:length(candrowidx))];
      % HKH: Update logical vector of candidates, candRow
      candRow(pivrow) = 0;
      candrowidx = find( candRow ); 
   end
end

% MODIFICATION LOG:
%
% 981112 mbk Rewritten for new input format. 
% 981112 hkh Some changes in comments
% 030108 hkh Change redundant constraints to -inf - inf constraints
% 030108 hkh Improve efficiency, remove bottlenecks, major revision
% 030109 hkh Add some handling of integer variables
% 030123 hkh Print summary if Prob.PriLevOpt > 0, more info if > 1
% 040125 hkh Use Prob.xxx(:) to safe guard for input row vectors
% 040306 hkh Correction of IntVars definition
% 040308 hkh Print Prob.b_L and b_L if PriLev > 1, (not > 0)
% 040422 hkh Return reduced problem
% 040325 hkh Set Prob.mLin if reduced problem
% 050330 ang Safeguard against empty ix in r_rw_sng
% 050330 hkh Put all statements with empty ix in r_rw_sng inside block 
% 070222 hkh Revised IntVars input format and definition of logical IV
