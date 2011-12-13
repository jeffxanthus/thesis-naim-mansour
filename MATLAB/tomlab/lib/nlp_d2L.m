% TOMLAB gateway routine
% Used as interface to filterSQP and MINLP
%
% nlp_d2L computes the Hessian
%
%       d2f(x) - lam' * d2c(x)
%
% to the Lagrangian function,
%
%   L(x,lam) =   f(x) - lam' * c(x)
%
% If input argument phase == 1, then only the second part lam'*d2c
% is calculated. If phase == 2, the full Hessian is calculated
%
% function d2L=nlp_d2L(x, lam, phase, Prob)
%
% Additionally, if giving a negative value for the phase parameter, the
% sign will be reversed in the Lagrangian. For example, if calling nlp_d2L 
% with d2L = nlp_d2L(x,lam,-2,Prob),  H(x) + lam'*d2c(x) will be returned.

% Anders Göran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.9.0$
% Written Nov 11, 2002.   Last modified Dec 22, 2009.

function d2L=nlp_d2L(x, lam, phs, Prob)

global n_H n_d2c NARG

global NLP_lamd2L NLP_xd2L NLP_d2L  % Global variables to store last d2L

sgn=sign(phs);
phase=abs(phs);

x=x(:);
d2H = [];

if ~isempty(NLP_xd2L) && ~isempty(NLP_lamd2L)
   if length(x)~=length(NLP_xd2L)
      NLP_xd2L=[];
   elseif length(lam)~=length(NLP_lamd2L)
      NLP_lamd2L=[];
   elseif all(x==NLP_xd2L) & all(lam==NLP_lamd2L)
      d2L=NLP_d2L;
      return
   end
end

if(isempty(lam))
   % If no nonlinear c/s
   d2L = nlp_H(x,Prob);
else
   
   if phase==1 % Only second part of Hessian
      d2L = nlp_d2c(x,lam,Prob);
   elseif phase==2 % Hessian of full Lagrangian
      d2H = nlp_H(x,Prob);
      if isempty(d2H)
         d2L = -sgn*nlp_d2c(x,lam,Prob);
      else
         d2L = d2H - sgn*nlp_d2c(x,lam,Prob);
      end
   end
end

d2L = sparse(d2L);

% if Prob.Warning & ~isempty(lam)
if Prob.Warning & Prob.probType ~= 2
   if(any(find(d2L-d2L')))
     % max(max(abs(d2L-d2L')))
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

NLP_xd2L = x;
NLP_lamd2L = lam;
NLP_d2L = d2L;

% MODIFICATION LOG:
%
% 021111 ango Wrote file
% 021112 ango Added handling of sparse Hessians
% 021223 hkh  Faster to avoid tests before doing sparse/full/transpose
% 021229 hkh  Change comments
% 021231 ango Always return sparse Hessian and transpose removed
% 030127 hkh  Display warning if not symmetric
% 040412 hkh  Only check if symmetric if Prob.Warning true
% 041016 hkh  Write long text and stop computations if unsymmetric d2L
% 041017 hkh  Correcting symmetry of d2L if MAD is used
% 050223 frhe x, lambda and d2L now saved in global variables
% 050902 ango Can do f+lam*c or f-lam*c - transparently
% 091222 hkh  Avoid test on symmetry if only H matrix (QP problem)
