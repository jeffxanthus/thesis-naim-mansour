% TOMLAB gateway routine
% Used as interface to BOEING Nonlinear Least Squares Solver
%
% nlp_LSd2L computes the Hessian
%
%       r' * d2r(x) - lam' * d2c(x)
%
% If input argument phase == 1, then only the second part lam'*d2c
% is calculated. If phase == 2, the full Hessian is calculated
%
% function LSd2L=nlp_LSd2L(x, lam, phase, Prob)

% Anders Göran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Jan 21, 2005.   Last modified Aug 13, 2009.

function d2L=nlp_LSd2L(x, lam, phase, Prob)

x=x(:);
d2H = [];

if(isempty(lam))
   % If no nonlinear c/s
   d2L = nlp_d2r(x, Prob);
else
   if phase==1 % Only second part of Hessian
     d2L = nlp_d2c(x,lam,Prob);
   elseif phase==2 % Hessian of full Lagrangian
     d2H = nlp_d2r(x, Prob);
     if isempty(d2H)
       d2L = - nlp_d2c(x,lam,Prob);
     else
       d2L = d2H - nlp_d2c(x,lam,Prob);
     end
   end
end

d2L = sparse(d2L);

if Prob.Warning
   if(any(find(d2L-d2L')))
    if Prob.ADCons == -1 | Prob.ADObj == -1
        d2L = 0.5*(d2L+d2L'); % ONLY IF MAD IS USED.
    else
      fprintf('\n\n')
      if isempty(lam)
         fprintf('The Hessian H of the objective is not symmetric!!!')
         fprintf('\n\n')
         fprintf('You must either correct your computations so that ')
         fprintf('the difference H-H'' is of the order of 1E-16 at most')
         fprintf('\n')
         fprintf('or make the Hessian symmetric by doing the following trick:')
         fprintf('\n')
         fprintf('H = 0.5*(H+H'');')
         fprintf('\n\n')
      elseif phase==1 | isempty(d2H)
         fprintf('The weighted Hessian d2L of the constraints ')
         fprintf('is not symmetric!!!')
         fprintf('\n\n')
         fprintf('You must either correct your computations so that ')
         fprintf('the difference d2L-d2L'' is of the order of 1E-16 at most')
         fprintf('\n')
         fprintf('or make d2L symmetric by doing the following trick:')
         fprintf('\n')
         fprintf('d2L = 0.5*(d2L+d2L'');')
         fprintf('\n\n')
      else
         if(any(find(d2H-d2H')))
            fprintf('The Hessian H of the objective is not symmetric!!!')
            fprintf('\n\n')
            fprintf('You must either correct your computations so that ')
            fprintf('the difference H-H'' is of the order of 1E-16 at most or')
            fprintf('\n')
            fprintf('make the Hessian symmetric by doing the following trick:')
            fprintf('\n')
            fprintf('H = 0.5*(H+H'');')
            fprintf('\n\n')
         end
         if(any(find((d2L-d2H)-(d2L'-d2H'))))
            fprintf('The weighted Hessian d2L of the constraints ')
            fprintf('is not symmetric!!!')
            fprintf('\n\n')
            fprintf('You must either correct your computations so that ')
            fprintf('the difference d2L-d2L'' is of the order of 1E-16 at most')
            fprintf('\n')
            fprintf('or make d2L symmetric by doing the following trick:')
            fprintf('\n')
            fprintf('d2L = 0.5*(d2L+d2L'');')
            fprintf('\n\n')
         end
      end
      error('Cannot proceed, Second order Lagrangian must be symmetric!!!');
    end
  end
end

% MODIFICATION LOG:
%
% 050121 frhe File created, based on nlp_d2L.m
% 090813 med  mlint check